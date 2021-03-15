!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz
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
! Add comments please!
!===================================================================================================================================

PROGRAM EstimateSigma_Bf
! MODULES
USE MOD_Globals
USE MOD_DG_Vars,                 ONLY: U
USE MOD_SwapMesh_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation,           ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Interpolation_Vars,      ONLY: xGP,wBary
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5,File_ID
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_SwapMesh,                ONLY: WriteNewStateFile,FinalizeSwapMesh
USE MOD_EOS,                     ONLY: ConsToPrim
USE MOD_InterpolateSolution,     ONLY: InterpolateSolution
USE MOD_MLMC_Vars
USE MOD_MLMC_SwapMesh_Vars
USE MOD_MLMC_Input,              ONLY: ReadSums,ReadSumsBF,ReadStateFile,ReadStateFileBF
USE MOD_MLMC_Output
USE MOD_IO_HDF5                 ,ONLY: AddToFieldData,InitMPIInfo
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,ReadAttribute
USE MOD_Mesh                    ,ONLY: DefineParametersMesh,InitMesh
USE MOD_Exactfunc               ,ONLY: ExactFunc
USE MOD_Analyze                 ,ONLY: InitAnalyzeBasis
USE MOD_EOS                     ,ONLY: DefineParametersEOS,InitEOS
USE MOD_Output_Vars             ,ONLY: ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)  :: tmp,H5PrmFile
INTEGER             :: i,j,k,iElem
REAL,ALLOCATABLE    :: UPrim(:,:,:,:,:)
LOGICAL             :: validInput,hasCoarse
INTEGER             :: nPrevious
!REAL                :: RefState(5)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
CALL InitMPIInfo()

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(1("*****************************************"))')
SWRITE(UNIT_stdOut,'(1("** START SIGMA ESTIMATE FOR BODYFORCES **"))')
SWRITE(UNIT_stdOut,'(1("*****************************************"))')
SWRITE(UNIT_stdOut,'(132("="))')

IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersInterpolation()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersEOS()
CALL DefineParametersEstSig()
CALL DefineParametersMesh()
CALL prms%CreateRealArrayOption("RefState"           , " RefState, for derived quantities")

! Parse parameters
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

! check if correct command line argumtns are given
validInput = .TRUE.
IF((nArgs.LT.4) .OR. (nArgs.GT.6))                   validInput=.FALSE.
IF(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))  validInput=.FALSE.
DO i=2,nArgs
  IF(.NOT.(STRICMP(GetFileExtension(Args(i)),'h5'))) validInput=.FALSE.
END DO
IF (.NOT.validInput) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: estimatesigma prm-file statefile_f bodyforces_f [statefile_c bodyforces_c]')
END IF

SWRITE(UNIT_stdOut,'(132("L"))')
CALL prms%read_options(Args(1))
SWRITE(UNIT_stdOut,'(132("L"))')

SWRITE(UNIT_stdOut,'(132("="))')

H5PrmFile=TRIM(Args(2))
StateFileFine=TRIM(Args(3))
BodyFocesFileFine=TRIM(Args(4))

CALL OpenDataFile(H5PrmFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'nPreviousRuns',1,IntScalar=nPrevious)
nStart=nPrevious+1
CALL ReadAttribute(File_ID,'nGlobalRuns',1,IntScalar=nEnd)
nEnd=nEnd+nPrevious
CALL CloseDataFile()

varAna       = GETINT('varAna')
! nAna         = GETINT('nAna')

hasCoarse = nArgs.EQ.6
IF(hasCoarse) THEN
  StateFileCoarse=TRIM(Args(5))
  BodyFocesFileCoarse=TRIM(Args(6))
END IF
!Swapmesh
! Read in polynomial degrees used for interpolation, supersampling and for the new state
! If they have not been specified, use the polynomial degree of the old state
NNew        = GETINT('NNew')
NInter      = NNew
NSuper      = NNew
useCurvedsOld = GETLOGICAL("useCurvedsOld",'.TRUE.')
useCurvedsNew = GETLOGICAL("useCurvedsNew",'.TRUE.')
printTroublemakers = GETLOGICAL('printTroublemakers','.TRUE.')
! Tolerance used to mark the position of an interpolation point in the new mesh as invalid. Will be invalid if the reference
! coordinate is more than maxtol outside of [-1,1], e.g "overshoot tolerance"
maxTol = 1. + GETREAL('maxTolerance','5.e-2')
! Set the abort tolerance, standard is the same as "overshoot tolerance"
WRITE(tmp,*) maxTol
abortTol     =GETREAL('abortTolerance',TRIM(tmp))
! New mesh file, the state will be interpolated to this one
MeshFileNew = GETSTR('MeshFileNew')

CALL InitInterpolation(NNew)
CALL InitIOHDF5()
#if USE_MPI
CALL InitMPIvars()
#endif

CALL InitSwapmesh(StateFileFine)
CALL SwapmeshToEstsig(.TRUE.)

IF (hasCoarse) THEN
  CALL FinalizeSwapMesh()
  CALL InitSwapmesh(StateFileCoarse)
  CALL SwapmeshToEstsig(.FALSE.)
END IF

CALL InitMesh(2,MeshFile_IN=MeshFileNew)
! CALL InitAnalyzeBasis(NNew,nAna,xGP,wBary)
CALL InitEos()

RefState = GETREALARRAY('RefState',nVar_State)

!-----------------------------------------------------------------------------------------------------------------------------------
CALL OpenDataFile(StateFileFine,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar=FileType)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL CloseDataFile()

IF (TRIM(FileType)=='State') THEN
  suffix='state'
  DataSetName="DG_Solution"
ELSEIF (TRIM(FileType)=='TimeAvg') THEN
  suffix='avg'
  DataSetName="Mean"
ELSE
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid FileType')
END IF
! TODO
DataSetNameBodyForces = 'BodyForces_BC_wall'

FileNameSums = 'postproc_'//TRIM(ProjectName)//'_'//TRIM(suffix)//'.h5'
FileNameSumsBF = 'postproc_'//TRIM(ProjectName)//'_'//'bodyforces'//'.h5'
snSamples  =1./REAL(nEnd)
snSamplesM1=1./(REAL(nEnd)-1.)

SWRITE(UNIT_stdOut,'(132("="))')

ALLOCATE(UFine        (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UFineSum     (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UFineSqSum   (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UCoarse      (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UCoarseSum   (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UCoarseSqSum (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(DUSqSum      (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(SigmaSqField (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(BiasField    (nVarTotal,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
IF (TRIM(FileType)=='State') THEN
  ALLOCATE(UPrim       (nVar_State+1,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
END IF

ALLOCATE(BodyForces             (1:9))
ALLOCATE(BodyForcesFine        (1:9))
ALLOCATE(BodyForcesFineSum     (1:9))
ALLOCATE(BodyForcesFineSqSum   (1:9))
ALLOCATE(BodyForcesCoarse      (1:9))
ALLOCATE(BodyForcesCoarseSum   (1:9))
ALLOCATE(BodyForcesCoarseSqSum (1:9))
ALLOCATE(DBodyForcesSqSum      (1:9))
ALLOCATE(SigmaSqBodyForces     (1:9))
ALLOCATE(BiasBodyForces        (1:9))

IF (.NOT.hasCoarse) UCoarse = 0.

UFineSum     =0.
UCoarseSum   =0.
UFineSqSum   =0.
UCoarseSqSum =0.
DUSqSum      =0.
BodyForcesFineSum = 0.
BodyForcesCoarseSum = 0.
BodyForcesFineSqSum = 0.
BodyForcesCoarseSqSum = 0.
DBodyForcesSqSum = 0.
nValBodyForces=(/9/)
IF (nStart.NE.1) THEN
  nVal=(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/)
  CALL ReadSums(FileNameSums)
  CALL ReadSumsBF(FileNameSumsBF)
END IF

!-----------------------------------------------------------------------------------------------------------------------------------

! Evaluate solution at new solution nodes
DO iSample=nStart,nEnd
  ! Read in the old state
  SWRITE(UNIT_stdOut,'(132("-"))')
  SWRITE(UNIT_stdOut,"(A,I0)") 'START PROCESSING SAMPLE ',iSample
  SWRITE(UNIT_stdOut,'(132("-"))')
  ! Read state file for each sample COARSE
  CALL ReadStateFileBF(BodyFocesFileFine,DataSetNameBodyForces,iSample-nStart)
  CALL EstsigToSwapmesh(.TRUE.)
  CALL ReadStateFile(StateFileFine,DataSetName,iSample-nStart)
  CALL InterpolateSolution()
  IF (TRIM(FileType)=='State') THEN
    DO iElem=1,nElemsNew
      DO k=0,NNew; DO j=0,NNew; DO i=0,NNew
        CALL ConsToPrim(UPrim(:,i,j,k,iElem),U(:,i,j,k,iElem))
      END DO; END DO; END DO
    END DO
    UFine(1:nVar_State,:,:,:,:)=U
    UFine(nVar_State+1:,:,:,:,:)=UPrim(2:,:,:,:,:)
    UFine(nVar_State*2+1,:,:,:,:)=(UPrim(5,:,:,:,:)-RefState(5))/(0.5*RefState(1)*NORM2(RefState(2:4))*NORM2(RefState(2:4)))
  ELSE
    UFine=U
  END IF
  BodyForcesFine = BodyForces
  IF (hasCoarse) THEN
    CALL ReadStateFileBF(BodyFocesFileCoarse,DataSetNameBodyForces,iSample-nStart)
    CALL EstsigToSwapmesh(.FALSE.)
    CALL ReadStateFile(StateFileCoarse,DataSetName,iSample-nStart)
    SWRITE(UNIT_stdOut,'(A)') ' EVALUATING THE COARSE SOLUTION ON NEW MESH ...'
    CALL InterpolateSolution()
    IF (TRIM(FileType)=='State') THEN
      DO iElem=1,nElemsNew
        DO k=0,NNew; DO j=0,NNew; DO i=0,NNew
          CALL ConsToPrim(UPrim(:,i,j,k,iElem),U(:,i,j,k,iElem))
        END DO; END DO; END DO
      END DO
      UCoarse(1:nVar_State,:,:,:,:)=U
      UCoarse(nVar_State+1:,:,:,:,:)=UPrim(2:,:,:,:,:)
      UCoarse(nVar_State*2+1,:,:,:,:)=(UPrim(5,:,:,:,:)-RefState(5))/(0.5*RefState(1)*NORM2(RefState(2:4))*NORM2(RefState(2:4)))
    ELSE
      UCoarse=U
    END IF
    BodyForcesCoarse = BodyForces
  END IF
  UFineSum     = UFineSum     + UFine
  UCoarseSum   = UCoarseSum   + UCoarse
  UFineSqSum   = UFineSqSum   + UFine*UFine
  UCoarseSqSum = UCoarseSqSum + UCoarse*UCoarse
  DUSqSum      = DUSqSum      + (UFine-UCoarse)*(UFine-UCoarse)

  BodyForcesFineSum     = BodyForcesFineSum     + BodyForcesFine
  BodyForcesCoarseSum   = BodyForcesCoarseSum   + BodyForcesCoarse
  BodyForcesFineSqSum   = BodyForcesFineSqSum   + BodyForcesFine*BodyForcesFine
  BodyForcesCoarseSqSum = BodyForcesCoarseSqSum + BodyForcesCoarse*BodyForcesCoarse
  DBodyForcesSqSum      = DBodyForcesSqSum      + (BodyForcesFine-BodyForcesCoarse)*(BodyForcesFine-BodyForcesCoarse)
END DO

SWRITE(UNIT_stdOut,'(132("="))')
! Bias:
BiasBodyForces = ABS( snSamples*(BodyForcesFineSum-BodyForcesCoarseSum))
Bias= BiasBodyForces(varAna)
! Stochastic Error:
SigmaSqBodyForces = snSamplesM1*( DBodyForcesSqSum - snSamples* (BodyForcesFineSum-BodyForcesCoarseSum) * (BodyForcesFineSum  -BodyForcesCoarseSum)  )
SigmaSq = SigmaSqBodyForces(varAna)
SigmaSqBodyForces = snSamplesM1*( BodyForcesFineSqSum - snSamples* (BodyForcesFineSum) * (BodyForcesFineSum )  )
SigmaSqFine = SigmaSqBodyForces(varAna)
CALL WriteSumsToHDF5()
CALL WriteBodyForcesSumsToHDF5()
CALL FinalizeSwapMesh()
CALL FinalizeEstimateSigma()

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' SigmaSq TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM EstimateSigma_Bf



!==================================================================================================================================
!> Define parameters for the used eos
!==================================================================================================================================
SUBROUTINE DefineParametersEstSig()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("EstimateSigmaSq")
CALL prms%CreateIntOption(      "varAna"             , "Specifiy the variable from U cons to comput SigmaSq")
CALL prms%CreateIntOption(      "NNew"               , "Polynomial degree used in new state files")
CALL prms%CreateIntOption(      "nAna"               , "Polynomial degree used for integration")
CALL prms%CreateStringOption(   "MeshFileOld"        , "Old mesh file (different than the one found in the state file)")
CALL prms%CreateStringOption(   "MeshFileNew"        , "New mesh file")
CALL prms%CreateLogicalOption(  "useCurvedsOld"      , "Controls usage of high-order information in old mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateLogicalOption(  "useCurvedsNew"      , "Controls usage of high-order information in new mesh. Turn off to discard "//&
                                                      "high-order data and treat curved meshes as linear meshes.", '.TRUE.')
CALL prms%CreateIntOption(      "NNew"               , "Polynomial degree used in new state files")
CALL prms%CreateRealOption(     "maxTolerance"       , "Tolerance used to mark points as invalid if outside of reference element "//&
                                                       "more than maxTolerance",'5.e-2')
CALL prms%CreateLogicalOption(  "printTroublemakers" , "Turn output of not-found points on or off",'.TRUE.')
CALL prms%CreateRealArrayOption("RefState"           , "If a RefState is defined, this state will be used at points that are "// &
                                                        "not found - without a RefState, the program will abort in this case")
CALL prms%CreateRealOption(     "abortTolerance"     , "Tolerance used to decide if the program should abort if no "// &
                                                       "RefState is given")
END SUBROUTINE DefineParametersEstSig



!===================================================================================================================================
!> Read in user defined parameters and prepare data for swapmesh.
!> The old and new mesh will be read and stored, the necessary Vandermonde matrizes are build and the parametric coordinates
!> of the new gauss points in the old mesh are found.
!===================================================================================================================================
SUBROUTINE InitSwapmesh(FileName)
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_SwapMesh_Vars
USE MOD_SwapMesh,                ONLY: ReadMeshCoords,PrepareVandermonde
USE MOD_ReadInTools
USE MOD_Commandline_Arguments
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension,INTTOSTR
USE MOD_DG_Vars,                 ONLY: U
USE MOD_Interpolation,           ONLY: GetVandermonde
USE MOD_Interpolation_Vars,      ONLY: NodeTypeVISU,NodeTypeCL
USE MOD_ChangeBasis,             ONLY: ChangeBasis3D
USE MOD_HDF5_Input,              ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute
USE MOD_IO_HDF5,                 ONLY: File_ID
USE MOD_SMParametricCoordinates, ONLY: GetParametricCoordinates
USE MOD_Interpolation,           ONLY: InitInterpolation
USE MOD_Mesh_Vars,               ONLY: nElems,OffsetElem,nGlobalElems
USE MOD_HDF5_Input,              ONLY: GetDataSize,DataSetExists
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN)  :: FileName       !< Mesh file on Level i to be read
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: i
INTEGER             :: nElems_State
REAL                :: Time
!===================================================================================================================================
! Open the first statefile to read necessary attributes
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=MeshFileOld)
CALL GetDataProps(nVar_State,NState,nElems_State,NodeTypeState)
CALL CloseDataFile()

! Initialize the old mesh, store the mesh coordinates (transformed to CL points) and the number of elements as well as the old NGeo
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(A)') ' INIT OLD MESH ...'
CALL ReadMeshCoords(MeshFileOld,useCurvedsOld,NGeoOld,nElemsOld,xCLOld)
SWRITE(UNIT_stdOut,*)'done in ',FLEXITIME()-Time

! Translate the old mesh along the displacement vector if needed
IF (CountOption('displacement').GE.1) THEN
  DO i=1,3
    xCLOld(i,:,:,:,:) = xCLOld(i,:,:,:,:)+displacement(i)
  END DO
END IF

! Initialize new mesh
Time=FLEXITIME()
SWRITE(UNIT_stdOut,'(A)') ' INIT NEW MESH ...'
CALL ReadMeshCoords(MeshFileNew,useCurvedsNew,NGeoNew,nElemsNew,xCLNew)
SWRITE(UNIT_stdOut,*)'done in ',FLEXITIME()-Time

! Set offset elem and local and global number of elements in mesh vars (later needed for output routine)
nGlobalElems = nElemsNew
nElems       = nElemsNew
OffsetElem   = 0 ! OffsetElem is 0 since the tool only works on singel

! Prepare the necessary Vandermonde matrizes, also interpolate the new mesh coordinates to xCLInter (on polynomial degree of NInter)
ALLOCATE(xCLInter(3,0:NInter,0:NInter,0:NInter,nElemsNew))
CALL prepareVandermonde()

! Evaluate parametric coordinates
SWRITE(UNIT_stdOut,'(A)') ' EVALUATING PARAMETRIC COORDINATES ...'
ALLOCATE(xiInter(3,  0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(InterToElem(0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(IPDone(     0:NInter,0:NInter,0:NInter,nElemsNew))
ALLOCATE(equalElem(nElemsNew))

! Find the parametric coordinates of the interpolation points of the new state in the old mesh
CALL GetParametricCoordinates()

ALLOCATE(U   (nVar_State,0:NNew,  0:NNew,  0:NNew,  nElemsNew))

END SUBROUTINE InitSwapmesh



!===================================================================================================================================
!> Copy all necessary variables for InterpolateSolution from InitSwapmesh variables to Estsig data structs
!===================================================================================================================================
SUBROUTINE SwapmeshToEstsig(isFine)
! MODULES
! !
USE MOD_SwapMesh_Vars
USE MOD_MLMC_Vars
USE MOD_MLMC_SwapMesh_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)      :: isFine
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(isFine) THEN
  FineNState              = NState
  FineNElemsOld           = NElemsOld

  ALLOCATE(FinexiInter(3,  0:NInter,0:NInter,0:NInter,nElemsNew))
  FinexiInter             = xiInter
  ALLOCATE(FineInterToElem(0:NInter,0:NInter,0:NInter,nElemsNew))
  FineInterToElem         = InterToElem
  ALLOCATE(FineIPDone(     0:NInter,0:NInter,0:NInter,nElemsNew))
  FineIPDone              = IPDone
  ALLOCATE(FineequalElem(nElemsNew))
  FineequalElem           = equalElem

  ALLOCATE(FineVdm_CLNInter_GPNNew(0:NNew,0:NInter))
  FineVdm_CLNInter_GPNNew = Vdm_CLNInter_GPNNew
  IF(NNew.NE.NState)THEN
    ALLOCATE(FineVdm_GPNState_GPNNew(0:NNew,0:NState))
    FineVdm_GPNState_GPNNew = Vdm_GPNState_GPNNew
  END IF
ELSE
  CoarseNState              = NState
  CoarseNElemsOld           = NElemsOld
  ALLOCATE(CoarsexiInter(3,  0:NInter,0:NInter,0:NInter,nElemsNew))
  CoarsexiInter             = xiInter
  ALLOCATE(CoarseInterToElem(0:NInter,0:NInter,0:NInter,nElemsNew))
  CoarseInterToElem         = InterToElem
  ALLOCATE(CoarseIPDone(     0:NInter,0:NInter,0:NInter,nElemsNew))
  CoarseIPDone              = IPDone
  ALLOCATE(CoarseequalElem(nElemsNew))
  CoarseequalElem           = equalElem

  ALLOCATE(CoarseVdm_CLNInter_GPNNew(0:NNew,0:NInter))
  CoarseVdm_CLNInter_GPNNew = Vdm_CLNInter_GPNNew
  IF(NNew.NE.NState)THEN
    ALLOCATE(CoarseVdm_GPNState_GPNNew(0:NNew,0:NState))
    CoarseVdm_GPNState_GPNNew = Vdm_GPNState_GPNNew
  END IF
END IF
END SUBROUTINE SwapmeshToEstsig



!===================================================================================================================================
!> Copy all necessary variables for InterpolateSolution from Estsig data type to according swapmesh readable variables
!===================================================================================================================================
SUBROUTINE EstsigToSwapmesh(isFine)
! MODULES

USE MOD_SwapMesh_Vars
USE MOD_MLMC_Vars
USE MOD_MLMC_SwapMesh_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
LOGICAL,INTENT(IN)      :: isFine
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
IF(isFine) THEN
  NState              = FineNState
  NElemsOld           = FineNElemsOld

  IF(NNew.NE.NState)THEN
    SDEALLOCATE(Vdm_GPNState_GPNNew)
    ALLOCATE(Vdm_GPNState_GPNNew(0:NNew,0:NState))
    Vdm_GPNState_GPNNew = FineVdm_GPNState_GPNNew
  END IF
  SDEALLOCATE(UOld)
  ALLOCATE(UOld(nVar_State,0:NState,0:NState,0:NState,nElemsOld))

  Vdm_CLNInter_GPNNew = FineVdm_CLNInter_GPNNew
  xiInter             = FinexiInter
  InterToElem         = FineInterToElem
  equalElem           = FineequalElem
  IPDone              = FineIPDone
ELSE
  NState              = CoarseNState
  NElemsOld           = CoarseNElemsOld

  IF(NNew.NE.NState)THEN
    SDEALLOCATE(Vdm_GPNState_GPNNew)
    ALLOCATE(Vdm_GPNState_GPNNew(0:NNew,0:NState))
    Vdm_GPNState_GPNNew = CoarseVdm_GPNState_GPNNew
  END IF
  SDEALLOCATE(UOld)
  ALLOCATE(UOld(nVar_State,0:NState,0:NState,0:NState,nElemsOld))

  Vdm_CLNInter_GPNNew = CoarseVdm_CLNInter_GPNNew
  xiInter             = CoarsexiInter
  InterToElem         = CoarseInterToElem
  equalElem           = CoarseequalElem
  IPDone              = CoarseIPDone
END IF
END SUBROUTINE EstsigToSwapmesh



SUBROUTINE IntegrateOverMesh(UIn,resu)
!===================================================================================================================================
! Integrate SigmaSq over mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: sJ
USE MOD_MLMC_Vars,          ONLY: nAna
USE MOD_Swapmesh_Vars,      ONLY: nElemsNew,NNew
USE MOD_Exactfunc,          ONLY: ExactFunc
USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Analyze_Vars,       ONLY: Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Vdm,FV_w
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)                 :: UIn(1,0:NNew,0:NNew,0:ZDIM(NNew),nElemsNew)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,k,l,m
REAL                            :: J_NAnalyze(1,0:nAna,0:nAna,0:ZDIM(nAna))
REAL                            :: U_NAnalyze(1,0:nAna,0:nAna,0:ZDIM(nAna), nElemsNew)
REAL                            :: J_N(1,0:NNew,0:NNew,0:ZDIM(NNew))
REAL                            :: IntegrationWeight
!===================================================================================================================================
resu=0.
DO iElem=1,nElemsNew
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,:,:,:)=1./sJ(:,:,:,iElem,0)
   CALL ChangeBasisVolume(1,NNew,nAna,Vdm_GaussN_NAnalyze,J_N(1:1,:,:,:),J_NAnalyze(1:1,:,:,:))
   CALL ChangeBasisVolume(1,NNew,nAna,Vdm_GaussN_NAnalyze,UIn(1:1,:,:,:,iElem),U_NAnalyze(1:1,:,:,:,iElem))
#if PP_dim ==3
   DO m=0,nAna
#else
   DO m=0,0
#endif
     DO l=0,nAna
       DO k=0,nAna
         IntegrationWeight = wGPVolAnalyze(k,l,m)*J_NAnalyze(1,k,l,m) !times 2 is a 2D-Fix, since sum over wGP_z=2
         resu = resu + U_NAnalyze(1,k,l,m,iElem)*IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO
END SUBROUTINE IntegrateOverMesh


SUBROUTINE BiasL2(UIn,resu)
!===================================================================================================================================
! Integrate SigmaSq over mesh
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: sJ
USE MOD_MLMC_Vars,          ONLY: nAna
USE MOD_Swapmesh_Vars,      ONLY: nElemsNew,NNew
USE MOD_Exactfunc,          ONLY: ExactFunc
USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Analyze_Vars,       ONLY: Vdm_GaussN_NAnalyze
USE MOD_Analyze_Vars,       ONLY: wGPVolAnalyze
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems,FV_Vdm,FV_w
#endif

! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)                 :: UIn(1,0:NNew,0:NNew,0:ZDIM(NNew),nElemsNew)
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(OUT)                :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: iElem,k,l,m
REAL                            :: J_NAnalyze(1,0:nAna,0:nAna,0:ZDIM(nAna))
REAL                            :: U_NAnalyze(1,0:nAna,0:nAna,0:ZDIM(nAna), nElemsNew)
REAL                            :: J_N(1,0:NNew,0:NNew,0:ZDIM(NNew))
REAL                            :: IntegrationWeight
!===================================================================================================================================
resu=0.
DO iElem=1,nElemsNew
   ! Interpolate the Jacobian to the analyze grid: be carefull we interpolate the inverse of the inverse of the jacobian ;-)
   J_N(1,:,:,:)=1./sJ(:,:,:,iElem,0)
   CALL ChangeBasisVolume(1,NNew,nAna,Vdm_GaussN_NAnalyze,J_N(1:1,:,:,:),J_NAnalyze(1:1,:,:,:))
   CALL ChangeBasisVolume(1,NNew,nAna,Vdm_GaussN_NAnalyze,UIn(1:1,:,:,:,iElem),U_NAnalyze(1:1,:,:,:,iElem))
   DO m=0,ZDIM(nAna)
     DO l=0,nAna
       DO k=0,nAna
         IntegrationWeight = wGPVolAnalyze(k,l,m)*J_NAnalyze(1,k,l,m) !times 2 is a 2D-Fix, since sum over wGP_z=2
         resu = resu + U_NAnalyze(1,k,l,m,iElem)*U_NAnalyze(1,k,l,m,iElem)*IntegrationWeight
       END DO ! k
     END DO ! l
   END DO ! m
END DO
resu = sqrt(resu)
END SUBROUTINE BiasL2


!===================================================================================================================================
!> Copy all necessary variables for InterpolateSolution from Estsig data type to according swapmesh readable variables
!===================================================================================================================================
SUBROUTINE FinalizeEstimateSigma()
! MODULES

USE MOD_MLMC_Vars
USE MOD_MLMC_SwapMesh_Vars
USE MOD_DG_Vars, ONLY: U,UPrim
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(U)
SDEALLOCATE(UPrim)

SDEALLOCATE(UFine)
SDEALLOCATE(UFineSum)
SDEALLOCATE(UFineSqSum)

SDEALLOCATE(UCoarse)
SDEALLOCATE(UCoarseSum)
SDEALLOCATE(UCoarseSqSum)

SDEALLOCATE(DUSqSum)
SDEALLOCATE(SigmaSqField)
!----------------------------------------------------------------------
SDEALLOCATE(FinexiInter)
SDEALLOCATE(FineInterToElem)
SDEALLOCATE(FineIPDone)
SDEALLOCATE(FineequalElem)
SDEALLOCATE(FineVdm_CLNInter_GPNNew)
SDEALLOCATE(FineVdm_GPNState_GPNNew)

SDEALLOCATE(CoarsexiInter)
SDEALLOCATE(CoarseInterToElem)
SDEALLOCATE(CoarseIPDone)
SDEALLOCATE(CoarseequalElem)
SDEALLOCATE(CoarseVdm_CLNInter_GPNNew)
SDEALLOCATE(CoarseVdm_GPNState_GPNNew)

END SUBROUTINE FinalizeEstimateSigma
