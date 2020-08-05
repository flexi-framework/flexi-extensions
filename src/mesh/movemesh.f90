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

!===================================================================================================================================
!> Contains subroutines to perform the movement of the the mesh.
!===================================================================================================================================
MODULE MOD_MoveMesh
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER :: MESHMOVETYPE_NONE            = 0
INTEGER,PARAMETER :: MESHMOVETYPE_CONSTANT_ACC    = 1
INTEGER,PARAMETER :: MESHMOVETYPE_CONSTANT        = 2
INTEGER,PARAMETER :: MESHMOVETYPE_ROTATION        = 3
INTEGER,PARAMETER :: MESHMOVETYPE_SLMI            = 4
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMoveMesh
  MODULE PROCEDURE InitMoveMesh
END INTERFACE

INTERFACE MoveMesh
  MODULE PROCEDURE MoveMesh
END INTERFACE

INTERFACE InterpolateMeshVel
  MODULE PROCEDURE InterpolateMeshVel
END INTERFACE

INTERFACE FinalizeMoveMesh
  MODULE PROCEDURE FinalizeMoveMesh
END INTERFACE

PUBLIC::InitMoveMesh,MoveMesh,InterpolateMeshVel,FinalizeMoveMesh,DefineParametersMoveMesh
!===================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters  of Mesh Movement
!==================================================================================================================================
SUBROUTINE DefineParametersMoveMesh()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("MoveMesh")
CALL prms%CreateLogicalOption(      'DoMeshAccelerate',"Use constant acceleration instead of constant velocity",multiple=.TRUE.)
CALL prms%CreateRealArrayOption(    'MeshVel',      "Mesh velocity for constant mesh movement",multiple=.TRUE.)
CALL prms%CreateRealArrayOption(    'MeshAcc',      "Constant mesh acceleration (DoMeshAccelerate = T)" ,multiple=.TRUE.)
CALL prms%CreateRealOption(         'RotationAngVel',"Angular velocity of rotation",multiple=.TRUE.)
END SUBROUTINE DefineParametersMoveMesh

!===================================================================================================================================
!> Get necessary parameters for mesh movement and initialize arrays.
!===================================================================================================================================
SUBROUTINE InitMoveMesh(MeshMode)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_ReadInTools
USE MOD_MoveMesh_Vars
USE MOD_MoveMesh_Mortar,      ONLY: InitMoveMeshMortar
USE MOD_Mesh_Vars,            ONLY: nSides,nElems,NGeo,MeshInitIsDone,NodeCoords
USE MOD_Mesh_Vars,            ONLY: IAmARotProc,IAmAStatProc,mySMPartition
USE MOD_Interpolation,        ONLY: GetVandermonde
USE MOD_Interpolation_Vars,   ONLY: NodeType,NodeTypeCL
USE MOD_IO_HDF5,              ONLY: AddToFieldData,FieldOut
USE MOD_Mappings,             ONLY: buildMappings
#if FV_ENABLED
USE MOD_Interpolation,        ONLY: GetNodesAndWeights
USE MOD_Basis,                ONLY: BarycentricWeights, InitializeVandermonde
USE MOD_FV_Basis,             ONLY: FV_Build_X_w_BdryX, FV_GetVandermonde,FV_Build_Vdm_Gauss_FVboundary
#endif
USE MOD_SM_Vars,              ONLY: SlidingMeshType,SlidingMeshDirection,nSlidingMeshPartitions
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN),OPTIONAL     :: MeshMode      !< Optional switch for POSTI
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=16)               :: DataSetName
CHARACTER(LEN=255)              :: VarNames(3)
INTEGER                         :: nVal(4)
#if FV_ENABLED
REAL                            :: xCL(0:NGeo), wBaryCL(0:NGeo)
REAL                            :: FV_w, FV_X(0:PP_N), FV_BdryX(0:PP_N+1)
#endif
INTEGER                         :: MeshModeLoc
INTEGER                         :: iPartition
REAL,PARAMETER                  :: eps=1.e-14
!===================================================================================================================================
MeshModeLoc = MERGE(MeshMode,1,PRESENT(MeshMode))

IF((.NOT.MeshInitIsDone).OR.MeshMoveInitIsDone)&
    CALL abort(__STAMP__,'InitMoveMesh not ready to be called or already called.')

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT MOVE MESH...'

ALLOCATE(DoMeshAccelerate(nSlidingMeshPartitions))
ALLOCATE(MeshVel(       3,nSlidingMeshPartitions))
ALLOCATE(MeshAcc(       3,nSlidingMeshPartitions))
ALLOCATE(RotationAngVel(  nSlidingMeshPartitions))
DoMeshAccelerate = .FALSE.
MeshAcc          = 0.
MeshVel          = 0.
RotationAngVel   = 0.

! Get mesh velocities
DO iPartition=1,nSlidingMeshPartitions
  SELECT CASE (SlidingMeshType(iPartition))
    CASE(SM_TYPE_PLANAR)
      ! We have to decide between motion with constant velocity and constant acceleration
      DoMeshAccelerate(iPartition) = GETLOGICAL('DoMeshAccelerate','F')

      IF(DoMeshAccelerate(iPartition)) THEN
        MeshAcc(:,iPartition)     = GETREALARRAY('MeshAcc',3,'0.,0.,0.')
        ! Currently only movement parallel to x/y-axis implemented
        IF((ABS(MeshAcc(MOD(SlidingMeshDirection(iPartition)*2,3),iPartition)+ABS(MeshVel(3,iPartition)))).GE.1e-15) THEN
          CALL abort(__STAMP__,'Input parameter MeshVel is not parallel to sliding mesh interface found in mesh file. & 
                               &The entries at following indices must be 0:'&
                               ,MOD(SlidingMeshDirection(iPartition)*2,3),REAL(3))
        END IF
        ! Set values manually to avoid parasitic velocities to best extent
        MeshAcc(MOD(SlidingMeshDirection(iPartition)*2,3),iPartition) = 0.
        MeshAcc(3                                        ,iPartition) = 0.
      ELSE
        MeshVel(:,iPartition)     = GETREALARRAY('MeshVel',3,'0.,0.,0.')
        ! Currently only movement parallel to x/y-axis implemented
        IF((ABS(MeshVel(MOD(SlidingMeshDirection(iPartition)*2,3),iPartition)+ABS(MeshVel(3,iPartition)))).GE.1e-15) THEN
          CALL abort(__STAMP__,'Input parameter MeshVel is not parallel to sliding mesh interface found in mesh file. & 
                               &The entries at following indices must be 0:'&
                               ,MOD(SlidingMeshDirection(iPartition)*2,3),REAL(3))
        END IF
        ! Set values manually to avoid parasitic velocities to best extent
        MeshVel(MOD(SlidingMeshDirection(iPartition)*2,3),iPartition) = 0.
        MeshVel(3                                        ,iPartition) = 0.
      END IF
    CASE(SM_TYPE_ROTATION_RADIAL,SM_TYPE_ROTATION_AXIAL)
      RotationAngVel(iPartition) = GETREAL('RotationAngVel')
    CASE DEFAULT
      CALL abort(__STAMP__,'Unknown SlidingMeshType in InitMoveMesh.')
  END SELECT
END DO

! Select correct MeshMoveType according to SlidingMeshType
IF(.NOT.IAmAStatProc) THEN ! Only Rot
  SELECT CASE(SlidingMeshType(mySMPartition))
    CASE(SM_TYPE_PLANAR)
      IF(DoMeshAccelerate(mySMPartition)) THEN
        MeshMoveType   = MESHMOVETYPE_CONSTANT_ACC
      ELSE
        MeshMoveType   = MESHMOVETYPE_CONSTANT
      END IF
    CASE(SM_TYPE_ROTATION_RADIAL,SM_TYPE_ROTATION_AXIAL)
      MeshMoveType   = MESHMOVETYPE_ROTATION
    CASE DEFAULT
      CALL abort(__STAMP__,'Unknown SlidingMeshType in InitMoveMesh.')
  END SELECT
ELSEIF(.NOT.IAmARotProc) THEN ! Only Stat 
  MeshMoveType = MESHMOVETYPE_NONE
ELSEIF(IAmAStatProc.AND.IAmARotProc) THEN ! Stat and Rot
#if USE_MPI
  CALL abort(__STAMP__,'With MPI a Proc should only have moving OR stationary elements.')
#else
  MeshMoveType = MESHMOVETYPE_SLMI
#endif
END IF

! Initialize actual Node position
ALLOCATE(NodeCoords_actual(3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems))
NodeCoords_actual=NodeCoords

! Initialize velocities on face and volume gauss/volume mesh points
ALLOCATE(v_NGeo(3,0:Ngeo,0:Ngeo,0:ZDIM(NGeo),nElems))
v_NGeo = 0.
ALLOCATE(v_NGeo_Face(1:3,0:NGeo,0:ZDIM(NGeo),1:nSides))
v_NGeo_Face = 0.
ALLOCATE(Elem_vGP(3,0:PP_N,0:PP_N,0:PP_NZ,nElems))
Elem_vGP = 0.
ALLOCATE(Face_vGP(3,0:PP_N,0:PP_NZ,1:nSides))
Face_vGP = 0.

! Initialize Vandermonde to interpolate mesh velcites from CL Ngeo to Gauss N
ALLOCATE(Vdm_CLNGeo_GaussN(0:PP_N,0:NGeo))
CALL GetVandermonde(Ngeo,NodeTypeCL,PP_N,NodeType,Vdm_CLNGeo_GaussN,modal=.FALSE.)

! Get mappings on NGeo
CALL buildMappings(NGeo,S2V2=S2V2_NGeo,FS2M=FS2M_NGeo,dim=PP_dim)

#if FV_ENABLED
ALLOCATE(Vdm_CLNGeo_FVBdryx(0:PP_N+1,0:NGeo))
ALLOCATE(Vdm_CLNGeo_FV     (0:PP_N  ,0:NGeo))
CALL GetNodesAndWeights(NGeo,NodeTypeCL,xCL)
CALL BarycentricWeights(NGeo, xCL, wBaryCL)
CALL FV_Build_X_w_BdryX(PP_N, FV_X, FV_w, FV_BdryX)
CALL InitializeVandermonde(NGeo,PP_N+1,wBaryCL,xCL,FV_BdryX,Vdm_CLNGeo_FVBdryx)
CALL InitializeVandermonde(NGeo,PP_N,  wBaryCL,xCL,FV_X,    Vdm_CLNGeo_FV)
#endif

! Add current mesh position to HDF5 output
DataSetName = 'Mesh_coordinates'
VarNames(1) = 'Mesh_coordinatesX'
VarNames(2) = 'Mesh_coordinatesY'
VarNames(3) = 'Mesh_coordinatesZ'
nVal        = (/3,NGeo+1,NGeo+1,ZDIM(NGeo)+1/)
CALL AddToFieldData(FieldOut,nVal,DataSetName,VarNames,NodeCoords_actual,doSeparateOutput=.TRUE.)

! Add current velocities to HDF5 output
DataSetName = 'Mesh_velocity'
VarNames(1) = 'Mesh_velocityX'
VarNames(2) = 'Mesh_velocityY'
VarNames(3) = 'Mesh_velocityZ'
nVal        = (/3,NGeo+1,NGeo+1,ZDIM(NGeo)+1/)
CALL AddToFieldData(FieldOut,nVal,DataSetName,VarNames,v_NGeo,doSeparateOutput=.TRUE.)

! Mortar specifics
!IF (MeshModeLoc.GT.0) CALL InitMoveMeshMortar()

MeshMoveInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT MOVE MESH DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitMoveMesh

!===================================================================================================================================
!> \brief Executes the movement of the mesh depending on the choosen type of movement.
!>
!> For each mesh points, the displacement and velocites are calculated based on the chosen mesh movement type.
!> Based on the new positions, the current metric terms will then be calculated. The Jacobian is either also calculated using the
!> mesh positions or will be updated by the GCL routines (needed to ensure free stream preservation on moving meshes).
!> The interpolation of the velocites  to the volume and surface gauss points will be done later in the DG routine since we need to
!> know if a DG/DG or a kind of FV interface is present.
!===================================================================================================================================
SUBROUTINE MoveMesh(t)
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_MoveMesh_Mortar,     ONLY : MortarMeshmovement
USE MOD_Interpolation_Vars,  ONLY : NodeType
USE MOD_MoveMesh_Vars,       ONLY : MeshMoveType,Elem_vGP,Face_vGP,MeshVel,MeshAcc,NodeCoords_actual
USE MOD_MoveMesh_Vars,       ONLY : RotationAngVel,RotationCenter,v_NGeo,MeshVelAlreadySet
USE MOD_MoveMesh_Vars,       ONLY : doMeshAccelerate
USE MOD_Mesh_Vars,           ONLY : RotatingElem,mySMPartition
USE MOD_Mesh_Vars,           ONLY : NGeo,nElems,NodeCoords
USE MOD_Mesh_Vars,           ONLY : Elem_xGP,Face_xGP,Elem_xGP_ref,Face_xGP_ref
USE MOD_Metrics,             ONLY : CalcMetrics,BuildCoords
USE MOD_ChangeBasisByDim,    ONLY : ChangeBasisVolume
#if FV_ENABLED
USE MOD_FV_Metrics,          ONLY : FV_CalcMetrics
USE MOD_FV_Vars,             ONLY : FV_Face_vGPXI,FV_Face_vGPEta
#if PP_dim == 3
USE MOD_FV_Vars,             ONLY : FV_Face_vGPZeta
#endif
#endif
USE MOD_SM_Vars,             ONLY : SlidingMeshDirection,SlidingMeshBoundaries
USE MOD_SM_Vars,             ONLY : SlidingMeshType
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)         :: t    !< Current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j,k,iElem,iPartition
REAL                    :: x(3)
REAL                    :: Displacement(1:3,0:NGeo,0:NGeo,0:ZDIM(NGeo),1:nElems) ! Temporary variable for displacement of each mesh
                                                                               ! node
REAL                    :: Period,rVec(2),phi,r
REAL                    :: myMeshVel,myMeshAcc
INTEGER                 :: myDirection
!===================================================================================================================================
! Select the type of mesh movement
SELECT CASE(MeshMoveType)
  CASE(MESHMOVETYPE_NONE) ! No movement
    ! ==============================================================================================================================
    ! TODO: Due to the current implementation of CalcMetrics, all metrics are communicated globally followed by a MPI_BARRIER.
    ! If any Sliding Mesh Partition needs to recalculate the metrics, all have to, otherwise MPI aborts. (Expensive!!)
    ! Improvement: introduce a local MPI communicator for each SM Partition. Until then we have to keep this ugly if-statement...
    IF (ANY((SlidingMeshType.EQ.SM_TYPE_ROTATION_RADIAL).OR.(SlidingMeshType.EQ.SM_TYPE_ROTATION_AXIAL))) THEN
      Displacement = 0.
    ELSE
      RETURN ! do nothing
    END IF

  CASE(MESHMOVETYPE_CONSTANT) ! constant velocity
    myDirection = SlidingMeshDirection(mySMPartition) ! Sliding Direction of current sm partition
    myMeshVel   = MeshVel(myDirection ,mySMPartition) ! MeshVel of current sm partition

    Period = (SlidingMeshBoundaries(2,mySMPartition)-SlidingMeshBoundaries(1,mySMPartition))

    ! calculate the new coordinates (they only change in sliding direction.)
    NodeCoords_actual(myDirection,:,:,:,:) = NodeCoords(  myDirection,:,:,:,:) + MOD(myMeshVel*t,Period) ! written to state file 
    Elem_xGP(         myDirection,:,:,:,:) = Elem_xGP_ref(myDirection,:,:,:,:) + MOD(myMeshVel*t,Period) ! interpolation points
    Face_xGP(         myDirection,:,:,:,:) = Face_xGP_ref(myDirection,:,:,:,:) + MOD(myMeshVel*t,Period) ! interpolation points face

    ! initialize the (constant) mesh velocities (only velocities in sliding direction)
    IF(.NOT.MeshVelAlreadySet) THEN
      v_NGeo(  myDirection,:,:,:,:) = myMeshVel ! written to state file for visu purposes
      Elem_vGP(myDirection,:,:,:,:) = myMeshVel ! interpolation points
      Face_vGP(myDirection,:,:,:)   = myMeshVel ! interpolation points face

#if FV_ENABLED    
      FV_Face_vGPXI  (myDirection,:,:,:,:) = myMeshVel
      FV_Face_vGPEta (myDirection,:,:,:,:) = myMeshVel
#if PP_dim==3       
      FV_Face_vGPZeta(myDirection,:,:,:,:) = myMeshVel
#endif
#endif /* FV_ENABLED */

      ! Set the velocities only once
      MeshVelAlreadySet = .TRUE.
    END IF

    ! ==============================================================================================================================
    ! TODO: Due to the current implementation of CalcMetrics, all metrics are communicated globally followed by a MPI_BARRIER.
    ! If any Sliding Mesh Partition needs to recalculate the metrics, all have to, otherwise MPI aborts. (Expensive!!)
    ! Improvement: introduce a local MPI communicator for each SM Partition. Until then we have to keep this ugly if-statement...
    IF (ANY((SlidingMeshType.EQ.SM_TYPE_ROTATION_RADIAL).OR.(SlidingMeshType.EQ.SM_TYPE_ROTATION_AXIAL))) THEN
      Displacement(          :,:,:,:,:) = 0.
      Displacement(myDirection,:,:,:,:) = MOD(myMeshVel*t,Period)
    ELSE
      ! Translatoric movement doesn't change the metrics
      RETURN
    END IF
  CASE(MESHMOVETYPE_CONSTANT_ACC) ! constant acceleration
    myDirection = SlidingMeshDirection(mySMPartition)       ! Sliding Direction of current sm partition
    myMeshAcc   = MeshAcc(myDirection ,mySMPartition)       ! MeshAcc of current sm partition

    Period = (SlidingMeshBoundaries(2,mySMPartition)-SlidingMeshBoundaries(1,mySMPartition))

    ! calculate the new coordinates (they only change in sliding direction.) s = 1/2*a*t^2
    NodeCoords_actual(myDirection,:,:,:,:) = NodeCoords(  myDirection,:,:,:,:) + MOD(0.5*myMeshAcc*t**2,Period) ! written to state file 
    Elem_xGP(         myDirection,:,:,:,:) = Elem_xGP_ref(myDirection,:,:,:,:) + MOD(0.5*myMeshAcc*t**2,Period) ! interpolation points
    Face_xGP(         myDirection,:,:,:,:) = Face_xGP_ref(myDirection,:,:,:,:) + MOD(0.5*myMeshAcc*t**2,Period) ! interpolation points face

    ! set the current mesh velocities (only in sliding direction) v = a*t
    v_NGeo(  myDirection,:,:,:,:) = myMeshAcc*t ! written to state file for visu purposes
    Elem_vGP(myDirection,:,:,:,:) = myMeshAcc*t ! interpolation points
    Face_vGP(myDirection,:,:,:)   = myMeshAcc*t ! interpolation points face*t*t

#if FV_ENABLED    
    FV_Face_vGPXI  (myDirection,:,:,:,:) = myMeshAcc*t
    FV_Face_vGPEta (myDirection,:,:,:,:) = myMeshAcc*t
#if PP_dim==3       
    FV_Face_vGPZeta(myDirection,:,:,:,:) = myMeshAcc*t
#endif
#endif /* FV_ENABLED */

    ! ==============================================================================================================================
    ! TODO: Due to the current implementation of CalcMetrics, all metrics are communicated globally followed by a MPI_BARRIER.
    ! If any Sliding Mesh Partition needs to recalculate the metrics, all have to, otherwise MPI aborts. (Expensive!!)
    ! Improvement: introduce a local MPI communicator for each SM Partition. Until then we have to keep this ugly if-statement...
    IF (ANY((SlidingMeshType.EQ.SM_TYPE_ROTATION_RADIAL).OR.(SlidingMeshType.EQ.SM_TYPE_ROTATION_AXIAL))) THEN
      Displacement(          :,:,:,:,:) = 0.
      Displacement(myDirection,:,:,:,:) = MOD(0.5*myMeshAcc*t**2,Period)
    ELSE
      ! Translatoric movement doesn't change the metrics
      RETURN
    END IF
  CASE(MESHMOVETYPE_ROTATION) ! Rigid rotation around z axis
    ! Period of rotation
    Period = 2.*PP_PI
    DO iElem=1,nElems
      DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
        x(:) = NodeCoords(:,i,j,k,iElem)
        ! Vector leading from rotation center to current point
        rVec(1) = x(1)-RotationCenter(1,mySMPartition)
        rVec(2) = x(2)-RotationCenter(2,mySMPartition)
        r = NORM2(rVec)
        ! Calculate angle
        phi = ATAN2(rVec(2),rVec(1))
        ! Add angle depending on current time
        phi = phi + MOD(t*RotationAngVel(mySMPartition),Period)
        ! Calculate new node coordinates
        Displacement(1,i,j,k,iElem) = RotationCenter(1,mySMPartition) + COS(phi)*r - x(1)
        Displacement(2,i,j,k,iElem) = RotationCenter(2,mySMPartition) + SIN(phi)*r - x(2)
        Displacement(3,i,j,k,iElem) = 0.
        ! Calculate mesh velocity
        v_NGeo(1,i,j,k,iElem) = -1.*SIN(phi)*r*RotationAngVel(mySMPartition)
        v_NGeo(2,i,j,k,iElem) =     COS(phi)*r*RotationAngVel(mySMPartition)
        v_NGeo(3,i,j,k,iElem) = 0.
      END DO; END DO; END DO! i,j,k=0,NGeo
    END DO ! iElem
   CASE(MESHMOVETYPE_SLMI) ! elements with different MESHMOVETYPE
     DO iElem=1,nElems
       ! Determine sm paritition of iElem
       iPartition = RotatingElem(iElem)
       SELECT CASE(SlidingMeshType(iPartition))
         CASE(SM_TYPE_NONE)
           Displacement(:,:,:,:,iElem) = 0.
           v_NGeo(      :,:,:,:,iElem) = 0.
         CASE(SM_TYPE_PLANAR)
           Period      = (SlidingMeshBoundaries(2,iPartition)-SlidingMeshBoundaries(1,iPartition))
           IF (DoMeshAccelerate(iPartition)) THEN ! constant acceleration
             DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
               ! Calculate new node coordinates
               Displacement(:,i,j,k,iElem) = MOD((t**2)*MeshAcc(:,iPartition),Period)
               ! Calculate mesh velocity
               v_NGeo(      :,i,j,k,iElem) = MeshAcc(:,iPartition)*t
             END DO; END DO; END DO! i,j,k=0,NGeo
           ELSE ! constant velocity
             DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
               ! Calculate new node coordinates
               Displacement(:,i,j,k,iElem) = MOD(t*MeshVel(:,iPartition),Period)
               ! Calculate mesh velocity
               v_NGeo(      :,i,j,k,iElem) = MeshVel(:,iPartition)
             END DO; END DO; END DO! i,j,k=0,NGeo
           END IF
         CASE(SM_TYPE_ROTATION_RADIAL,SM_TYPE_ROTATION_AXIAL)
           ! Period of rotation
           Period = 2.*PP_PI
           DO k=0,ZDIM(NGeo); DO j=0,NGeo; DO i=0,NGeo
             x(:) = NodeCoords(:,i,j,k,iElem)
             ! Vector leading from rotation center to current point
             rVec(1) = x(1)-RotationCenter(1,iPartition)
             rVec(2) = x(2)-RotationCenter(2,iPartition)
             r = NORM2(rVec)
             ! Calculate angle
             phi = ATAN2(rVec(2),rVec(1))
             ! Add angle depending on current time
             phi = phi + MOD(t*RotationAngVel(iPartition),Period)
             ! Calculate new node coordinates
             Displacement(1,i,j,k,iElem) = RotationCenter(1,iPartition) + COS(phi)*r - x(1)
             Displacement(2,i,j,k,iElem) = RotationCenter(2,iPartition) + SIN(phi)*r - x(2)
             Displacement(3,i,j,k,iElem) = 0.
             ! Calculate mesh velocity
             v_NGeo(1,i,j,k,iElem) = -1.*SIN(phi)*r*RotationAngVel(iPartition)
             v_NGeo(2,i,j,k,iElem) =     COS(phi)*r*RotationAngVel(iPartition)
             v_NGeo(3,i,j,k,iElem) = 0.
           END DO; END DO; END DO! i,j,k=0,NGeo
         CASE DEFAULT
           CALL Abort(__STAMP__, &
                    'Unknown SlidingMeshType in MoveMesh.')
       END SELECT
     END DO ! iElem
  CASE DEFAULT
      CALL Abort(__STAMP__, &
               'This case is not supported in movemesh.')
END SELECT

! For mortar meshes: Get a continuous representation of the displacement and velocites on the faces
!CALL MortarMeshmovement(Displacement,v_NGeo)

! Add displacement to node positions -> current node coordinates
NodeCoords_actual = NodeCoords + Displacement

! Calc new metrics and also new Elem_xGP
#if GCL
CALL CalcMetrics(NodeCoords_actual(:,:,:,:,:),InitJacobian=.FALSE.)
#else
CALL CalcMetrics(NodeCoords_actual(:,:,:,:,:),InitJacobian=.TRUE.)
#endif /*GCL*/
CALL BuildCoords(NodeCoords_actual,NodeType,PP_N,Elem_xGP)

#if FV_ENABLED
CALL FV_CalcMetrics()  ! FV metrics
#endif

END SUBROUTINE MoveMesh

!===================================================================================================================================
!> \brief Interpolate the mesh velocities from mesh representation to volume and surface gauss points.
!> 
!> In this routine, the mesh point velocites computed by the MoveMesh routines are interpolated to the volume and surface gauss
!> points, where they are needed to compute the fluxes. Depending if the cell is a DG or a FV cell, the volume operation is done
!> by a simple ChangeBasis call or by using a routine that converts the velocities to the FV representation.
!===================================================================================================================================
SUBROUTINE InterpolateMeshVel(v_CLNGeo,v_CLNGeo_Face)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: NGeo,nElems,nSides
USE MOD_MoveMesh_Vars    ,ONLY: Elem_vGP,Vdm_CLNGeo_GaussN,MeshMoveType
#if FV_ENABLED
USE MOD_FV_Vars          ,ONLY: FV_Elems
#endif
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_ChangeBasisByDim, ONLY: ChangeBasisVolume
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)    :: v_CLNGeo(     3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems) !< volume mesh velocity
REAL,INTENT(INOUT)    :: v_CLNGeo_Face(3,0:NGeo,0:ZDIM(NGeo),nSides)        !< face   mesh velocity
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iElem
!===================================================================================================================================
IF ((MeshMoveType.LE.MESHMOVETYPE_CONSTANT).AND.(.NOT. postiMode)) RETURN ! For no or constant mesh movement, no interpolation is necessary

DO iElem=1,nElems
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN ! DG element
#endif
    ! Get mesh velocities on every Gauss point by change of basis - CL_NGeo to Gauss N
    CALL ChangeBasisVolume(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,v_CLNGeo(:,:,:,:,iElem),Elem_vGP(:,:,:,:,iElem))
#if FV_ENABLED
  ELSE ! FV element
    ! Interpolate to FV representation
    CALL InterpolateToFV(v_CLNGeo(:,:,:,:,iElem),iElem)
  END IF
#endif
END DO
! Interpolate to face values
CALL InterpolateToFaces(v_CLNGeo,v_CLNGeo_Face)
END SUBROUTINE InterpolateMeshVel

#if FV_ENABLED
!===================================================================================================================================
!> Interpolate mesh velocity from cheb.-lob. NGeo volume data calculated by move mesh routine to volume data in Finite Volume
!> representation.
!===================================================================================================================================
SUBROUTINE InterpolateToFV(v_NGeo,iElem)
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: NGeo
USE MOD_MoveMesh_Vars    ,ONLY: Vdm_CLNGeo_FV,Vdm_CLNGeo_FVBdryx
USE MOD_ChangeBasis      ,ONLY: ChangeBasis1D
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_FV_Vars          ,ONLY: FV_Face_vGPXI,FV_Face_vGPEta
#if PP_dim == 3
USE MOD_FV_Vars          ,ONLY: FV_Face_vGPZeta
#endif
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER,INTENT(IN) :: iElem                               !< local element ID
REAL,INTENT(INOUT) :: v_NGeo(3,0:NGeo,0:NGeo,0:ZDIM(NGeo))  !< volume mesh velocity of the element
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: i,j,k
REAL               :: tmp(3,0:NGeo,0:ZDIM(NGeo),0:PP_N+1)
!===================================================================================================================================
! Xi direction
DO k=0,ZDIM(NGeo); DO j=0,NGeo
  CALL ChangeBasis1D  (3,NGeo,PP_N+1,Vdm_CLNGeo_FVBdryx,v_NGeo(:,:,j,k),tmp(:,j,k,:))
END DO; END DO! j,k
DO i=1,PP_N
  CALL ChangeBasisSurf(3,NGeo,PP_N,Vdm_CLNGeo_FV,tmp(:,:,:,i),FV_Face_vGPXI(:,:,:,i,iElem))
END DO

! Eta direction
DO k=0,ZDIM(NGeo); DO i=0,NGeo
  CALL ChangeBasis1D  (3,NGeo,PP_N+1,Vdm_CLNGeo_FVBdryx,v_NGeo(:,i,:,k),tmp(:,i,k,:))
END DO; END DO! i,k
DO j=1,PP_N
  CALL ChangeBasisSurf(3,NGeo,PP_N,Vdm_CLNGeo_FV,tmp(:,:,:,j),FV_Face_vGPEta(:,:,:,j,iElem))
END DO

#if (PP_dim == 3)
! Zeta direction
DO j=0,ZDIM(NGeo); DO i=0,NGeo
  CALL ChangeBasis1D  (3,NGeo,PP_N+1,Vdm_CLNGeo_FVBdryx,v_NGeo(:,i,j,:),tmp(:,i,j,:))
END DO; END DO! i,j
DO k=1,PP_N
  CALL ChangeBasisSurf(3,NGeo,PP_N,Vdm_CLNGeo_FV,tmp(:,:,:,k),FV_Face_vGPZeta(:,:,:,k,iElem))
END DO
#endif /* PP_dim == 3 */

END SUBROUTINE InterpolateToFV
#endif /* FV_ENABLED */


!===================================================================================================================================
!> Interpolate mesh velocity from cheb.-lob. NGeo volume data calculated by move mesh routine to Gauss N face data needed to
!> fill the fluxes.
!===================================================================================================================================
SUBROUTINE InterpolateToFaces(v_CLNGeo,v_CLNGeo_Face)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: ElemToSide,NGeo,SideToElem,nElems,nSides,lastMPISide_MINE
USE MOD_MoveMesh_Vars    ,ONLY: Face_vGP,Vdm_CLNGeo_GaussN,S2V2_NGeo
#if FV_ENABLED
USE MOD_MoveMesh_Vars    ,ONLY: Vdm_CLNGeo_FV
USE MOD_FV_Vars          ,ONLY: FV_Elems_sum
#endif
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(IN)    :: v_CLNGeo(     3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems) !< volume mesh velocity
REAL,INTENT(INOUT) :: v_CLNGeo_Face(3,0:NGeo,0:ZDIM(NGeo),nSides)        !< face   mesh velocity
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q,iLocSide,pq(2)
INTEGER :: SideID,flip
INTEGER :: iElem
REAL    :: tmp(3,0:NGeo,0:ZDIM(NGeo))
!===================================================================================================================================

! Loop over all sides that are attached to an element (mortar sides with no elements attached are already done in the mortar
! routine) and get the velocity on the face.
DO iElem=1,nElems
#if PP_dim == 3
  DO iLocSide=1,6
#else 
  DO iLocSide=2,5
#endif
    SideID=ElemToSide(E2S_SIDE_ID,iLocSide,iElem)
    ! Since the velocities are continuous across faces, they only need to be computed once per side. This is done by the master element
    ! if both elements belonging to the side are on this proc or by the slave otherwise (MPI or mortar).

         !      We are not master                  and          there is a master
    IF ((ElemToSide(E2S_FLIP,iLocSide,iElem).NE.0).AND.(SideToElem(S2E_ELEM_ID,SideID).GT.0)) CYCLE

    ! Select slice of volume data depending on the local side
    SELECT CASE(iLocSide)
    CASE(XI_MINUS)
      tmp=v_CLNGeo(1:3,0   ,:   ,:   ,iElem)
    CASE(XI_PLUS)
      tmp=v_CLNGeo(1:3,NGeo,:   ,:   ,iElem)
    CASE(ETA_MINUS)
      tmp=v_CLNGeo(1:3,:   ,0   ,:   ,iElem)
    CASE(ETA_PLUS)
      tmp=v_CLNGeo(1:3,:   ,NGeo,:   ,iElem)
    CASE(ZETA_MINUS)
      tmp=v_CLNGeo(1:3,:   ,:   ,0   ,iElem)
    CASE(ZETA_PLUS)
      tmp=v_CLNGeo(1:3,:   ,:   ,NGeo,iElem)
    END SELECT

    ! turn into right hand system of side
    flip = ElemToSide(E2S_FLIP,iLocSide,iElem)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,iLocSide)
      v_CLNGeo_Face(1:3,p,q,SideID)=tmp(:,pq(1),pq(2))
    END DO; END DO ! p,q
  END DO ! iLocSide=1,6

END DO ! iElem


! Loop over all sides we need to calculate the fluxes on. There the mesh velocity is needed!
! Convert the velocity based on NGeo to the solution gauss points.
DO SideID = 1, lastMPISide_MINE
!  IF ((SideID.GE.firstSMSide).AND.(SideID.LE.lastSMSide)) CYCLE
#if FV_ENABLED
  IF (FV_Elems_sum(SideID).EQ.0) THEN ! DG/DG interface
#endif /* FV_ENABLED */
    ! Change basis to Gauss points
    CALL ChangeBasisSurf(3,NGeo,PP_N,Vdm_CLNGeo_GaussN,v_CLNGeo_Face(:,:,:,SideID),Face_vGP(:,:,:,SideID))
#if FV_ENABLED
  ELSE ! FV/DG,DG/FV,FV/FV interface
    ! Convert Ngeo to FV, flux will be calculated on the FV subcell points 
    CALL ChangeBasisSurf(3,NGeo,PP_N,Vdm_CLNGeo_FV    ,v_CLNGeo_Face(:,:,:,SideID),Face_vGP(:,:,:,SideID))
  END IF
#endif
END DO ! iSide = 1, lastMPISide_MINE
END SUBROUTINE InterpolateToFaces

!===================================================================================================================================
!> Finalizes global variables used by MoveMesh
!===================================================================================================================================
SUBROUTINE FinalizeMoveMesh()
! MODULES                                                                                                                          !
USE MOD_MoveMesh_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
MeshMoveInitIsDone = .FALSE.
SDEALLOCATE(RotationCenter)
SDEALLOCATE(RotationAngVel)
SDEALLOCATE(MeshVel)
SDEALLOCATE(MeshAcc)
SDEALLOCATE(DoMeshAccelerate)
SDEALLOCATE(v_NGeo)
SDEALLOCATE(v_NGeo_Face)
SDEALLOCATE(Elem_vGP)
SDEALLOCATE(Face_vGP)
SDEALLOCATE(NodeCoords_actual)
SDEALLOCATE(Vdm_CLNGeo_GaussN)
SDEALLOCATE(S2V2_NGeo)
#if FV_ENABLED
SDEALLOCATE(Vdm_CLNGeo_FVBdryx)
SDEALLOCATE(Vdm_CLNGeo_FV)
#endif
SDEALLOCATE(SmallInnerMortarSideIDs)
SDEALLOCATE(NbProcsMortarSend)
SDEALLOCATE(nMortarSidesSend)
SDEALLOCATE(offsetMortarSidesSend)
SDEALLOCATE(mapMortarSidesSend)
SDEALLOCATE(SmallMortarMPISidesRcv)
SDEALLOCATE(NbProcsMortarRcv)
SDEALLOCATE(nMortarSidesRcv)
SDEALLOCATE(offsetMortarSidesRcv)
SDEALLOCATE(mapMortarSidesRcv)

END SUBROUTINE FinalizeMoveMesh

END MODULE MOD_MoveMesh
