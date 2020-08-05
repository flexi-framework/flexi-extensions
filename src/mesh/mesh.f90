!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

!==================================================================================================================================
!> Contains control routines to read high-order meshes, provide mesh data to the solver, build the metrics, partition the domain.
!==================================================================================================================================
MODULE MOD_Mesh
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES (PUBLIC)
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitMesh
  MODULE PROCEDURE InitMesh
END INTERFACE

INTERFACE FinalizeMesh
  MODULE PROCEDURE FinalizeMesh
END INTERFACE

PUBLIC::InitMesh
PUBLIC::FinalizeMesh
!==================================================================================================================================

PUBLIC::DefineParametersMesh
CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Mesh")
CALL prms%CreateStringOption(  'MeshFile',            "(relative) path to meshfile (mandatory).")
CALL prms%CreateLogicalOption( 'useCurveds',          "Controls usage of high-order information in mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption( 'interpolateFromTree', "For non-conforming meshes, built by refinement from a tree structure, "//&
                                                      "the metrics can be built from the tree geometry if it is contained "//&
                                                      "in the mesh. Can improve free-stream preservation.",&
                                                      '.TRUE.')
CALL prms%CreateRealOption(    'meshScale',           "Scale the mesh by this factor (shrink/enlarge).",&
                                                      '1.0')
CALL prms%CreateLogicalOption( 'meshdeform',          "Apply simple sine-shaped deformation on cartesion mesh (for testing).",&
                                                      '.FALSE.')
#if (PP_dim == 3)
CALL prms%CreateLogicalOption( 'crossProductMetrics', "Compute mesh metrics using cross product form. Caution: in this case "//&
                                                      "free-stream preservation is only guaranteed for N=3*NGeo.",&
                                                      '.FALSE.')
#endif
CALL prms%CreateIntOption(     'debugmesh',           "Output file with visualization and debug information for the mesh. "//&
                                                      "0: no visualization, 3: Paraview binary",'0')
CALL prms%CreateStringOption(  'BoundaryName',        "Names of boundary conditions to be set (must be present in the mesh!)."//&
                                                      "For each BoundaryName a BoundaryType needs to be specified.",&
                                                      multiple=.TRUE.)
CALL prms%CreateIntArrayOption('BoundaryType',        "Type of boundary conditions to be set. Format: (BC_TYPE,BC_STATE)",&
                                                      multiple=.TRUE.)
CALL prms%CreateLogicalOption( 'writePartitionInfo',  "Write information about MPI partitions into a file.",'.FALSE.')
CALL prms%CreateIntOption(     'NGeoOverride',        "Override switch for NGeo. Interpolate mesh to different NGeo." //&
                                                      "<1: off, >0: Interpolate",'-1')
END SUBROUTINE DefineParametersMesh


!==================================================================================================================================
!> Routine controlling the initialization of the mesh.
!> - parameter and mesh reading
!> - domain partitioning
!> - allocation of mesh arrays
!> - build mesh mappings to handle volume/surface operations
!> - compute the mesh metrics
!> - provide mesh metrics for overintegration
!==================================================================================================================================
SUBROUTINE InitMesh(meshMode,MeshFile_IN)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars
USE MOD_HDF5_Input
USE MOD_Interpolation_Vars, ONLY:InterpolationInitIsDone,NodeTypeVISU,NodeTypeCL,NodeType
USE MOD_ChangeBasisByDim   ,ONLY:ChangeBasisVolume
USE MOD_Interpolation      ,ONLY:GetVandermonde
USE MOD_Mesh_ReadIn,        ONLY:readMesh
USE MOD_Prepare_Mesh,       ONLY:setLocalSideIDs,fillMeshInfo
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETSTR,GETREAL,GETINT
USE MOD_Metrics,            ONLY:BuildCoords,CalcMetrics
USE MOD_DebugMesh,          ONLY:writeDebugMesh
USE MOD_Mappings,           ONLY:buildMappings
USE MOD_ChangeBasis,        ONLY:ChangeBasis3D
USE MOD_Interpolation,      ONLY:GetVandermonde
#if USE_MPI
USE MOD_Prepare_Mesh,       ONLY:exchangeFlip
#endif
#if FV_ENABLED
USE MOD_FV_Metrics,         ONLY:InitFV_Metrics,FV_CalcMetrics
#endif
USE MOD_IO_HDF5,            ONLY:AddToElemData,ElementOut
#if (PP_dim == 2)
USE MOD_2D
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: meshMode !< 0: only read and build Elem_xGP,
                               !< 1: as 0 + build connectivity, 2: as 1 + calc metrics
CHARACTER(LEN=255),INTENT(IN),OPTIONAL :: MeshFile_IN !< file name of mesh to be read
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL              :: x(3),meshScale
REAL,POINTER      :: coords(:,:,:,:,:),coordsTmp(:,:,:,:,:),Vdm_CLNGeo_CLNGeoOverride(:,:)
INTEGER           :: iElem,i,j,k,nElemsLoc
LOGICAL           :: validMesh
REAL,ALLOCATABLE  :: Vdm_EQNgeo_CLNgeo(:,:) ! Vandermonde from equidistant mesh points to CL points on NGeo
REAL,ALLOCATABLE  :: XCL_Ngeo(:,:,:,:)      ! CL mesh points in a single element
INTEGER           :: firstMasterSide     ! lower side ID of array U_master/gradUx_master...
INTEGER           :: lastMasterSide      ! upper side ID of array U_master/gradUx_master...
INTEGER           :: firstSlaveSide      ! lower side ID of array U_slave/gradUx_slave...
INTEGER           :: lastSlaveSide       ! upper side ID of array U_slave/gradUx_slave...
INTEGER           :: iSide,LocSideID,SideID
INTEGER           :: NGeoOverride
INTEGER           :: nLocMortars
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.MeshInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitMesh not ready to be called or already called.')
END IF

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A,I1,A)') ' INIT MESH IN MODE ',meshMode,'...'

! prepare pointer structure (get nElems, etc.)
IF (PRESENT(MeshFile_IN)) THEN
  MeshFile = MeshFile_IN
ELSE
  MeshFile = GETSTR('MeshFile')
END IF
validMesh = ISVALIDMESHFILE(MeshFile)
IF(.NOT.validMesh) &
    CALL CollectiveStop(__STAMP__,'ERROR - Mesh file not a valid HDF5 mesh.')

useCurveds=GETLOGICAL('useCurveds','.TRUE.')
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'Ngeo',1,IntScalar=NGeo)
CALL CloseDataFile()

IF(useCurveds.AND.(PP_N.LT.NGeo))THEN
  SWRITE(UNIT_stdOut,'(A)') 'WARNING: N<NGeo, for curved hexa normals are only approximated,&
                           & can cause problems on periodic boundaries! Set N>=NGeo'
ENDIF

CALL readMesh(MeshFile) !set nElems

! Preparation step: Interpolate mesh points form equidistant to CL on Ngeo.
! Subsequent operations from mesh movement will live on CL points!
ALLOCATE(Vdm_EQNgeo_CLNgeo(0:Ngeo,0:Ngeo))
CALL GetVandermonde(Ngeo,NodeTypeVISU,Ngeo,NodeTypeCL,Vdm_EQNgeo_CLNgeo,modal=.FALSE.)

ALLOCATE(XCL_Ngeo(3,0:Ngeo,0:Ngeo,0:Ngeo))
DO iElem=1,nElems
  CALL ChangeBasis3D(3,NGeo,NGeo,Vdm_EQNGeo_CLNGeo,NodeCoords(:,:,:,:,iElem),XCL_Ngeo)
  ! Overwrite NodeCoords array (equidistant up to now)
  NodeCoords(:,:,:,:,iElem) = XCL_Ngeo
END DO

DEALLOCATE(Vdm_EQNgeo_CLNgeo)
DEALLOCATE(XCL_Ngeo)

#if (PP_dim == 2)
! If this is a two dimensional calculation, all subsequent operations are performed on the reduced mesh.
SWRITE(UNIT_StdOut,'(A)') " RUNNING A 2D SIMULATION! "
! The mesh coordinates read in by the readMesh routine are therefore reduced by one dimension.
CALL to2D_rank5((/1,0,0,0,1/),(/3,NGeo,NGeo,NGeo,nElems/),4,NodeCoords)
NodeCoords(3,:,:,:,:) = 0.
#endif

! if trees are available: compute metrics on tree level and interpolate to elements
interpolateFromTree=.FALSE.
IF(isMortarMesh) interpolateFromTree=GETLOGICAL('interpolateFromTree','.TRUE.')

! At the moment, interpolation from tree not working for moving meshes!
IF (interpolateFromTree) &
      CALL abort(__STAMP__,'At the moment, interpolation from tree not working for moving meshes!')

IF(interpolateFromTree)THEN
#if (PP_dim == 2)
  CALL CollectiveStop(__STAMP__,&
      "interpolateFromTree not supported in 2D.")
#endif
  coords=>TreeCoords
  NGeo=NGeoTree
  nElemsLoc=nTrees
ELSE
  coords=>NodeCoords
  nElemsLoc=nElems
ENDIF

NGeoOverride=GETINT('NGeoOverride','-1')
IF(NGeoOverride.GT.0)THEN
  ALLOCATE(CoordsTmp(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
  ALLOCATE(Vdm_CLNGeo_CLNGeoOverride(0:NGeoOverride,0:NGeo))
  CALL GetVandermonde(Ngeo, NodeTypeCL, NgeoOverride, NodeTypeCL, Vdm_CLNgeo_CLNgeoOverride, modal=.FALSE.)
  DO iElem=1,nElemsLoc
    CALL ChangeBasisVolume(3,Ngeo,NgeoOverride,Vdm_CLNGeo_CLNGeoOverride,coords(:,:,:,:,iElem),coordsTmp(:,:,:,:,iElem))
  END DO
  ! cleanup
  IF(interpolateFromTree)THEN
    DEALLOCATE(TreeCoords)
    ALLOCATE(TreeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>TreeCoords
  ELSE
    DEALLOCATE(NodeCoords)
    ALLOCATE(NodeCoords(3,0:NGeoOverride,0:NGeoOverride,0:NGeoOverride,nElemsLoc))
    coords=>NodeCoords
  END IF
  Coords = CoordsTmp
  NGeo = NGeoOverride
  DEALLOCATE(CoordsTmp, Vdm_CLNGeo_CLNGeoOverride)
END IF

SWRITE(UNIT_StdOut,'(a3,a30,a3,i0)')' | ','Ngeo',' | ', Ngeo

! scale and deform mesh if desired (warning: no mesh output!)
meshScale=GETREAL('meshScale','1.0')
IF(ABS(meshScale-1.).GT.1e-14)&
  Coords = Coords*meshScale

IF(GETLOGICAL('meshdeform','.FALSE.'))THEN
  DO iElem=1,nElemsLoc
    DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
      x(:)=Coords(:,i,j,k,iElem)
#if PP_dim==3
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))*SIN(PP_Pi*x(3))
#else
      Coords(:,i,j,k,iElem) = x+ 0.1*SIN(PP_Pi*x(1))*SIN(PP_Pi*x(2))
#endif
    END DO; END DO; END DO;
  END DO
END IF

! Build the coordinates of the solution gauss points in the volume
ALLOCATE(Elem_xGP(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
IF(interpolateFromTree)THEN
  CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP,TreeCoords)
ELSE
  CALL BuildCoords(NodeCoords,NodeType,PP_N,Elem_xGP)
ENDIF

! Return if no connectivity and metrics are required (e.g. for visualization mode)
IF (meshMode.GT.0) THEN

  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING setLocalSideIDs..."
  CALL setLocalSideIDs()

#if USE_MPI
  ! for MPI, we need to exchange flips, so that MINE MPISides have flip>0, YOUR MpiSides flip=0
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING exchangeFlip..."
  CALL exchangeFlip()
#endif

  !RANGES
  !-----------------|-----------------|-------------------|
  !    U_master     | U_slave         |    FLUX           |
  !-----------------|-----------------|-------------------|
  !  BCsides        |                 |    BCSides        |
  !  InnerMortars   |                 |    InnerMortars   |
  !  InnerSides     | InnerSides      |    InnerSides     |
  !  MPI_MINE sides | MPI_MINE sides  |    MPI_MINE sides |
  !                 | MPI_YOUR sides  |    MPI_YOUR sides |
  !  MPIMortars     |                 |    MPIMortars     |
  !-----------------|-----------------|-------------------|

  firstBCSide          = 1
  firstSMSide          = firstBCSide         +nBCSides
  firstMortarInnerSide = firstSMSide         +nSMSides
  firstInnerSide       = firstMortarInnerSide+nMortarInnerSides
  firstMPISide_MINE    = firstInnerSide      +nInnerSides
  firstMPISide_YOUR    = firstMPISide_MINE   +nMPISides_MINE
  firstMortarMPISide   = firstMPISide_YOUR   +nMPISides_YOUR

  lastBCSide           = firstSMSide-1
  lastSMSide           = firstMortarInnerSide-1 
  lastMortarInnerSide  = firstInnerSide    -1
  lastInnerSide        = firstMPISide_MINE -1
  lastMPISide_MINE     = firstMPISide_YOUR -1
  lastMPISide_YOUR     = firstMortarMPISide-1
  lastMortarMPISide    = nSides


  firstMasterSide = 1
  lastMasterSide  = nSides
  firstSlaveSide  = firstInnerSide
  lastSlaveSide   = lastMPISide_YOUR
  nSidesMaster    = lastMasterSide-firstMasterSide+1
  nSidesSlave     = lastSlaveSide -firstSlaveSide+1

  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMasterSide     ', firstMasterSide,lastMasterSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastSlaveSide      ', firstSlaveSide, lastSlaveSide
  LOGWRITE(*,*)'-------------------------------------------------------'
  LOGWRITE(*,'(A25,I8,I8)')'first/lastBCSide         ', firstBCSide         ,lastBCSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastSMSide         ', firstSMSide         ,lastSMSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMortarInnerSide', firstMortarInnerSide,lastMortarInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastInnerSide      ', firstInnerSide      ,lastInnerSide
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_MINE   ', firstMPISide_MINE   ,lastMPISide_MINE
  LOGWRITE(*,'(A25,I8,I8)')'first/lastMPISide_YOUR   ', firstMPISide_YOUR   ,lastMPISide_YOUR
  LOGWRITE(*,'(A30,I8,I8)')'first/lastMortarMPISide  ', firstMortarMPISide  ,lastMortarMPISide
  LOGWRITE(*,*)'-------------------------------------------------------'

  ! fill ElemToSide, SideToElem,BC
  ALLOCATE(ElemToSide(2,6,nElems))
  ALLOCATE(SideToElem(5,nSides))
  ALLOCATE(BC(1:nBCSides))
  ALLOCATE(AnalyzeSide(1:nSides))
  ElemToSide  = 0
  SideToElem  = -1   !mapping side to elem, sorted by side ID (for surfint)
  BC          = 0
  AnalyzeSide = 0

  !NOTE: nMortarSides=nMortarInnerSides+nMortarMPISides
  ALLOCATE(MortarType(2,1:nSides))              ! 1: Type, 2: Index in MortarInfo
  ALLOCATE(MortarInfo(MI_FLIP,4,nMortarSides)) ! [1]: 1: Neighbour sides, 2: Flip, [2]: small sides
  MortarType=0
  MortarInfo=-1

  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING fillMeshInfo..."
  CALL fillMeshInfo()

  ! dealloacte pointers
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING deleteMeshPointer..."
  CALL deleteMeshPointer()

#if (PP_dim ==2)
  ! In 2D, there is only one flip for the slave sides (1)
  SideToElem(S2E_FLIP,:) = MIN(1,SideToElem(S2E_FLIP,:))
  ElemToSide(:,1,:) = -999
  ElemToSide(:,6,:) = -999
  ElemToSide(E2S_FLIP,:,:) = MIN(1,ElemToSide(E2S_FLIP,:,:))
  MortarInfo(MI_FLIP,:,:)  = MIN(1,MortarInfo(MI_FLIP,:,:))
#endif

  ! Build necessary mappings
  CALL buildMappings(PP_N,V2S=V2S,S2V=S2V,S2V2=S2V2,FS2M=FS2M,dim=PP_dim)
END IF

IF (meshMode.GT.1) THEN

  ! ----- CONNECTIVITY IS NOW COMPLETE AT THIS POINT -----

  ! volume data
  ALLOCATE(      dXCL_N(3,3,0:PP_N,0:PP_N,0:PP_NZ,nElems)) ! temp
  ALLOCATE(Metrics_fTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
  ALLOCATE(Metrics_gTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
  ALLOCATE(Metrics_hTilde(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
  ALLOCATE(            sJ(  0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
  ALLOCATE(     scaledJac(  0:PP_N,0:PP_N,0:PP_NZ,nElems))
  NGeoRef=3*NGeo ! build jacobian at higher degree
  ALLOCATE(    DetJac_Ref(1,0:NgeoRef,0:NgeoRef,0:ZDIM(NGeoRef),nElems))

  ! surface data
  ALLOCATE(      Face_xGP(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  ALLOCATE(       NormVec(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  ALLOCATE(      TangVec1(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  ALLOCATE(      TangVec2(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  ALLOCATE(      SurfElem(  0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  ALLOCATE(     Ja_Face(3,3,0:PP_N,0:PP_NZ,             1:nSides)) ! temp


  ! assign all metrics Metrics_fTilde,Metrics_gTilde,Metrics_hTilde
  ! assign 1/detJ (sJ)
  ! assign normal and tangential vectors and surfElems on faces

#if (PP_dim == 3)
  ! compute metrics using cross product instead of curl form (warning: no free stream preservation!)
  crossProductMetrics=GETLOGICAL('crossProductMetrics','.FALSE.')
#endif
  SWRITE(UNIT_stdOut,'(A)') "NOW CALLING calcMetrics..."
  CALL CalcMetrics(NodeCoords,InitJacobian=.TRUE.)     ! DG metrics
  
  ! Save reference GP position for better performance in moving mesh
  ALLOCATE(Elem_xGP_ref(3,0:PP_N,0:PP_N ,0:PP_NZ     ,1:nElems))
  ALLOCATE(Face_xGP_ref(3,0:PP_N,0:PP_NZ,0:FV_ENABLED,1:nSides))
  Elem_xGP_ref = Elem_xGP
  Face_xGP_ref = Face_xGP

#if FV_ENABLED
  CALL InitFV_Metrics()  ! FV metrics
  SWRITE(UNIT_stdOut,'(A)',ADVANCE='NO') '  Build FV-Metrics ...'
  CALL FV_CalcMetrics()  ! FV metrics
  SWRITE(UNIT_stdOut,'(A)')' Done !'
#endif
  ! debugmesh: param specifies format to output, 0: no output, 1: tecplot ascii, 2: tecplot binary, 3: paraview binary
  CALL WriteDebugMesh(GETINT('debugmesh','0'))
END IF

IF (meshMode.GT.0) THEN
  ALLOCATE(SideToGlobalSide(nSides))
  DO iElem=1,nElems
    nLocMortars=0

#if PP_dim == 3
    DO LocSideID=1,6
#else
    DO LocSideID=2,5
#endif
      SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
      iSide  = ElemInfo(3,iElem+offsetElem) + LocSideID + nLocMortars

      ! Skip side if it is a virtual small mortar (104: linear, 204: curved) and increment element mortar counter
      DO WHILE ((SideInfo(1,iSide).EQ.104).OR.(SideInfo(1,iSide).EQ.204))
        nLocMortars = nLocMortars + 1
        iSide  = ElemInfo(3,iElem+offsetElem) + LocSideID + nLocMortars
      END DO

      SideToGlobalSide(SideID) = ABS(SideInfo(2,iSide))
    END DO
  END DO ! iElem
END IF

SDEALLOCATE(TreeCoords)
SDEALLOCATE(xiMinMax)
SDEALLOCATE(ElemToTree)
! Do not deallocate this here as we may actually re-calcualted it for deforming meshes.
!IF (.NOT.postiMode) DEALLOCATE(scaledJac)

CALL AddToElemData(ElementOut,'myRank',IntScalar=myRank)

! Build SM communication structure 
IF(DoSlidingMesh) CALL SlidingMeshConnect(MeshMode,MeshFile,meshScale)

MeshInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitMesh



!==================================================================================================================================
!> This routine reads in the SlidingMeshInfo. In axial direction of the rotation, the sides on the sliding mesh interface are 
!> sorted per layer. In each layer, the nAzimuthalSides are sorted by an increasing angle. 
!> The SlidingMeshInfo gives the globalSideID for the n-th side in the i-th layer
!==================================================================================================================================
SUBROUTINE SlidingMeshConnect(MeshMode,FileString,meshScale)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_HDF5_Input,    ONLY:OpenDataFile,CloseDataFile,DatasetExists,File_ID,GetDataSize,HSize,nDims,ReadArray,ReadAttribute
USE MOD_Mesh_Vars,     ONLY:nAzimuthalSides,nLayers,SlidingMeshInfo,SideToGlobalSide,firstSMSide,lastSMSide,SideToElem &
                            ,AL_ToRotProc,AL_ToStatProc,IntToPart,IAmAStatProc,IAmARotProc,mySMPartition
USE MOD_MoveMesh_Vars, ONLY:RotationCenter
USE MOD_SM_Vars,       ONLY:SlidingMeshDirection,SlidingMeshBoundaries,nSlidingMeshInterfaces,nSlidingMeshPartitions
USE MOD_Mesh_Vars,     ONLY:nElems,RotatingElem,offsetElem
USE MOD_Mesh_Vars,     ONLY:nStatProcs,GlobalSlidingMeshInfo
USE MOD_Mesh_Vars,     ONLY:LocToGlobInterface,GlobToLocInterface
USE MOD_IO_HDF5,       ONLY:AddToElemData,ElementOut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)              :: meshScale  !< (IN) meshScale
INTEGER,INTENT(IN)           :: MeshMode   !< (IN) mesh mode flag
CHARACTER(LEN=*),INTENT(IN)  :: FileString !< (IN) mesh filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                      :: SlidingMeshInfoFound
INTEGER                      :: Offset=0 ! Every process reads all BCs
INTEGER                      :: iLayer,iAzimuthalSide
!INTEGER,ALLOCATABLE          :: GlobalSlidingMeshInfo(:,:)
INTEGER                      :: nGlobalSlidingMeshInterfaces
INTEGER                      :: GlobalSideID,iSide
INTEGER                      :: iInterface,iGlobInterface
CHARACTER(LEN=255)           :: str
INTEGER                      :: idx
#if USE_MPI
INTEGER                      :: iAz,iL
#endif
!==================================================================================================================================
! read-in Sliding Mesh Data
CALL OpenDataFile(FileString,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

CALL ReadAttribute(File_ID,'nSlidingMeshInterfaces',1,IntScalar=nGlobalSlidingMeshInterfaces)
CALL ReadAttribute(File_ID,'nSlidingMeshPartitions',1,IntScalar=nSlidingMeshPartitions)

ALLOCATE(GlobalSlidingMeshInfo(nGlobalSlidingMeshInterfaces))
ALLOCATE(nAzimuthalSides(nGlobalSlidingMeshInterfaces))
ALLOCATE(nLayers(        nGlobalSlidingMeshInterfaces))
ALLOCATE(IntToPart(      nGlobalSlidingMeshInterfaces))

ALLOCATE(SlidingMeshDirection(   1:nSlidingMeshPartitions))
ALLOCATE(RotationCenter(3,       1:nSlidingMeshPartitions))
ALLOCATE(SlidingMeshBoundaries(2,1:nSlidingMeshPartitions))

CALL ReadAttribute(File_ID,'SM_InterfaceToPartition',nGlobalSlidingMeshInterfaces,IntArray=IntToPart)

CALL ReadArray('SlidingMeshDirection', 1,(/  nSlidingMeshPartitions/),0,1, IntArray=SlidingMeshDirection)
CALL ReadArray('SlidingMeshCenter',    2,(/3,nSlidingMeshPartitions/),0,1,RealArray=RotationCenter)
CALL ReadArray('SlidingMeshBoundaries',2,(/2,nSlidingMeshPartitions/),0,1,RealArray=SlidingMeshBoundaries)
! Instantly multiply Boundaries with mesh scaling factor
SlidingMeshBoundaries = meshScale * SlidingMeshBoundaries

! read-in SlidingMeshInfo
DO iInterface=1,nGlobalSlidingMeshInterfaces
  WRITE(str,'(I0)') iInterface
  CALL DatasetExists(File_ID,TRIM('SlidingMeshInfo'//str),SlidingMeshInfoFound)
  IF(.NOT.SlidingMeshInfoFound) CALL CollectiveStop(&
    __STAMP__,&
        "SlidingMeshInfo not found in meshfile!")
  CALL GetDataSize(File_ID,TRIM('SlidingMeshInfo'//str),nDims,HSize)
  CHECKSAFEINT(HSize(1),4)
  nAzimuthalSides(iInterface)=INT(HSize(1),4)
  CHECKSAFEINT(HSize(2),4)
  nLayers(        iInterface)=INT(HSize(2),4)
  ALLOCATE(GlobalSlidingMeshInfo(iInterface)%Sides(1:nAzimuthalSides(iInterface),1:nLayers(iInterface)))
  OffSet=0
  CALL ReadArray(TRIM('SlidingMeshInfo'//str),2,(/nAzimuthalSides(iInterface),nLayers(iInterface)/),0,1,IntArray=GlobalSlidingMeshInfo(iInterface)%Sides)
END DO

!Debug output
ALLOCATE(RotatingElem(1:nElems))
CALL ReadArray('RotatingElem',1,(/nElems/),offsetElem,1,IntArray=RotatingElem)
CALL CloseDataFile()

! Only build SM mappings for meshMode greater then 0
IF (MeshMode.GT.0) THEN
  ! map the global-SideID from the SlidingMeshInfo to the local side structure
  ALLOCATE(SlidingMeshInfo(1:5,firstSMSide:lastSMSide))
  SlidingMeshInfo=0

  ! Get number of different SM Interfaces adjacent to elems of this proc
  ALLOCATE(GlobToLocInterface(1:nGlobalSlidingMeshInterfaces))
  GlobToLocInterface = 0
  nSlidingMeshInterfaces = 0
  DO iSide=firstSMSide,lastSMSide
    GlobalSideID=SideToGlobalSide(iSide)
    DO iInterface=1,nGlobalSlidingMeshInterfaces
      DO iLayer=1,nLayers(iInterface)
        DO iAzimuthalSide=1,nAzimuthalSides(iInterface)
          IF(GlobalSideID.EQ.GlobalSlidingMeshInfo(iInterface)%Sides(iAzimuthalSide,iLayer))THEN
            GlobToLocInterface(iInterface) = 1
          END IF
        END DO
      END DO
    END DO
  END DO
  nSlidingMeshInterfaces = SUM(GlobToLocInterface)

  ALLOCATE(LocToGlobInterface(1:nSlidingMeshInterfaces))
  LocToGlobInterface = 0
  ! Build mapping from local to global interface
  idx = 1
  DO iInterface=1,nGlobalSlidingMeshInterfaces
    IF(GlobToLocInterface(iInterface).EQ.1) THEN
      LocToGlobInterface(idx) = iInterface
      idx = idx + 1
    END IF
  END DO

  ! match each local sliding mesh side with it's azimuthal position and position in the i-th layer
  ! loop over all SM sides to map the side to it's position, layer, ElemID and Neighbor-ElemID
  DO iInterface=1,nSlidingMeshInterfaces
    iGlobInterface=LocToGlobInterface(iInterface)
    ! loop over all layers and azimuthal-sides per layer to get the position and layer of the side
    DO iSide=firstSMSide,lastSMSide
      GlobalSideID=SideToGlobalSide(iSide)
      DO iLayer=1,nLayers(iGlobInterface)
        DO iAzimuthalSide=1,nAzimuthalSides(iGlobInterface)
          IF(GlobalSideID.EQ.GlobalSlidingMeshInfo(iGlobInterface)%Sides(iAzimuthalSide,iLayer))THEN
            SlidingMeshInfo(SM_SIDE_TO_INTERFACE ,iSide)=iGlobInterface
            SlidingMeshInfo(SM_SIDE_TO_POSITION  ,iSide)=iAzimuthalSide
            SlidingMeshInfo(SM_SIDE_TO_LAYER     ,iSide)=iLayer
            ! and store the corresponding elements, duplicate but simpler?
            SlidingMeshInfo(SM_SIDE_TO_ELEM_ID   ,iSide)=SideToElem(S2E_ELEM_ID,iSide)
            SlidingMeshInfo(SM_SIDE_TO_NB_ELEM_ID,iSide)=SideToElem(S2E_NB_ELEM_ID,iSide)
          END IF
        END DO ! iAzimuthalSide=1,nAzimuthalSides
      END DO ! iLayer=1,nLayers
    END DO ! iSide=firstSMSide,lastSMSide
  END DO ! iInterface1,nSlidingMeshInterface
END IF

CALL AddToElemData(ElementOut,'RotatingElem',IntArray=RotatingElem)

! Deallocate data structure
DO iInterface=1,nGlobalSlidingMeshInterfaces
  SDEALLOCATE(GlobalSlidingMeshInfo(iInterface)%Sides)
END DO
SDEALLOCATE(GlobalSlidingMeshInfo)

! Save sm partition if proc only has moving elems
IF(IAmAStatProc) THEN
  mySMPartition = -1
ELSEIF(IAmARotProc) THEN
  mySMPartition = RotatingElem(1)
END IF

! TODO: Replace with data structure to avoid unneccessary big array
ALLOCATE(AL_ToStatProc(1:nGlobalSlidingMeshInterfaces,1:MAXVAL(nAzimuthalSides),1:MAXVAL(nLayers)))
ALLOCATE(AL_ToRotProc( 1:nGlobalSlidingMeshInterfaces,1:MAXVAL(nAzimuthalSides),1:MAXVAL(nLayers)))
AL_ToStatProc=-1
AL_ToRotProc=-1

!write a matrix with all sm sides sorted by azimuth and layer, 
!and communicate to all procs in which proc the side is
#if USE_MPI
DO iSide=firstSMSide,lastSMSide
  DO iInterface=1,nSlidingMeshInterfaces
    iGlobInterface = LocToGlobInterface(iInterface)
    IF (iGlobInterface.NE.SlidingMeshInfo(SM_SIDE_TO_INTERFACE  ,iSide)) CYCLE
    iAz=SlidingMeshInfo(SM_SIDE_TO_POSITION  ,iSide)
    iL =SlidingMeshInfo(SM_SIDE_TO_LAYER     ,iSide)
    IF(myRank.LT.nStatProcs)THEN
      AL_ToStatProc(iGlobInterface,iAz,iL)=myRank
    ELSE
      AL_ToRotProc( iGlobInterface,iAz,iL)=myRank
    END IF
  END DO
END DO
CALL MPI_ALLREDUCE(MPI_IN_PLACE,AL_ToStatProc,nGlobalSlidingMeshInterfaces*MAXVAL(nAzimuthalSides)*MAXVAL(nLayers),MPI_INTEGER,MPI_MAX,MPI_COMM_FLEXI,iError)
CALL MPI_ALLREDUCE(MPI_IN_PLACE,AL_ToRotProc, nGlobalSlidingMeshInterfaces*MAXVAL(nAzimuthalSides)*MAXVAL(nLayers),MPI_INTEGER,MPI_MAX,MPI_COMM_FLEXI,iError)
#else
AL_ToStatProc=0
AL_ToRotProc= 0
#endif
END SUBROUTINE SlidingMeshConnect



!============================================================================================================================
!> Deallocate mesh data.
!============================================================================================================================
SUBROUTINE FinalizeMesh()
! MODULES
USE MOD_Mesh_Vars
USE MOD_Mappings  ,ONLY:FinalizeMappings
#if FV_ENABLED
USE MOD_FV_Vars   ,ONLY:FV_Elems_master
USE MOD_FV_Metrics,ONLY:FinalizeFV_Metrics
#endif
IMPLICIT NONE
!============================================================================================================================
! Deallocate global variables, needs to go somewhere else later
SDEALLOCATE(ElemToSide)
SDEALLOCATE(SideToElem)
SDEALLOCATE(BC)
SDEALLOCATE(dXCL_N)
SDEALLOCATE(Ja_Face)
SDEALLOCATE(AnalyzeSide)

SDEALLOCATE(MortarType)
SDEALLOCATE(MortarInfo)


! allocated during ReadMesh
SDEALLOCATE(NodeCoords)
SDEALLOCATE(BoundaryName)
SDEALLOCATE(BoundaryType)


! Volume
SDEALLOCATE(Elem_xGP)
SDEALLOCATE(Elem_xGP_ref)
SDEALLOCATE(Metrics_fTilde)
SDEALLOCATE(Metrics_gTilde)
SDEALLOCATE(Metrics_hTilde)
SDEALLOCATE(sJ)
SDEALLOCATE(scaledJac)
SDEALLOCATE(DetJac_Ref)

! surface
SDEALLOCATE(Face_xGP)
SDEALLOCATE(Face_xGP_ref)
SDEALLOCATE(NormVec)
SDEALLOCATE(TangVec1)
SDEALLOCATE(TangVec2)
SDEALLOCATE(SurfElem)

! ijk sorted mesh
SDEALLOCATE(Elem_IJK)
SDEALLOCATE(ElemInfo)
SDEALLOCATE(SideInfo)
SDEALLOCATE(SideToGlobalSide)

! mappings
CALL FinalizeMappings()

!> sliding mesh
SDEALLOCATE(IntToPart)
SDEALLOCATE(nAzimuthalSides)
SDEALLOCATE(nLayers)
SDEALLOCATE(LocToGlobInterface)
SDEALLOCATE(GlobToLocInterface)
SDEALLOCATE(SlidingMeshInfo)
SDEALLOCATE(RotatingElem)
SDEALLOCATE(AL_ToStatProc)
SDEALLOCATE(AL_ToRotProc )
SDEALLOCATE(nRotElems)
SDEALLOCATE(nRotProcs)

#if FV_ENABLED
SDEALLOCATE(FV_Elems_master) ! moved here from fv.f90
CALL FinalizeFV_Metrics()
#endif

MeshInitIsDone = .FALSE.
END SUBROUTINE FinalizeMesh

END MODULE MOD_Mesh
