!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
#include "eos.h"

!==================================================================================================================================
!> This module contains the necessary routines for communication with Relexi via a SmartRedis client.
!==================================================================================================================================
MODULE MOD_SmartRedis

#if USE_SMARTREDIS
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER      :: PRM_SMARTREDIS_NONE      = 0
INTEGER,PARAMETER      :: PRM_SMARTREDIS_HIT       = 1
INTEGER,PARAMETER      :: PRM_SMARTREDIS_CHANNEL   = 2
INTEGER,PARAMETER      :: PRM_SMARTREDIS_CYLINDER  = 3

INTERFACE DefineParametersSmartRedis
  MODULE PROCEDURE DefineParametersSmartRedis
END INTERFACE

INTERFACE InitSmartRedis
  MODULE PROCEDURE InitSmartRedis
END INTERFACE

INTERFACE AnalyzeSmartRedis
  MODULE PROCEDURE AnalyzeSmartRedis
END INTERFACE

INTERFACE FinalizeSmartRedis
  MODULE PROCEDURE FinalizeSmartRedis
END INTERFACE

!INTERFACE GatheredWriteSmartRedis
!  MODULE PROCEDURE GatheredWriteSmartRedis
!END INTERFACE

INTERFACE ExchangeDataSmartRedis
  MODULE PROCEDURE ExchangeDataSmartRedis
END INTERFACE


PUBLIC :: DefineParametersSmartRedis
PUBLIC :: InitSmartRedis
PUBLIC :: FinalizeSmartRedis
PUBLIC :: ExchangeDataSmartRedis
PUBLIC :: AnalyzeSmartRedis
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for SmartRedis
!==================================================================================================================================
SUBROUTINE DefineParametersSmartRedis()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("SmartRedis")
CALL prms%CreateLogicalOption("doSmartRedis", "Communicate via the SmartRedis Client", ".FALSE.")
CALL prms%CreateIntFromStringOption('SR_Type', "Type of SmartRedis communication: None, HIT, Channel, Cylinder.", "hit")
CALL addStrListEntry('SR_Type','none',     PRM_SMARTREDIS_NONE)
CALL addStrListEntry('SR_Type','hit',      PRM_SMARTREDIS_HIT)
CALL addStrListEntry('SR_Type','channel',  PRM_SMARTREDIS_CHANNEL)
CALL addStrListEntry('SR_Type','cylinder', PRM_SMARTREDIS_CYLINDER)
CALL prms%CreateLogicalOption("SR_useInvariants", "Use Invariants of gradient tensor as state for agent", ".FALSE.")
CALL prms%CreateLogicalOption("SR_doNormInvariants", "Normalizing invariants of velocity gradient tensor", ".TRUE.")
CALL prms%CreateIntOption("SR_nVarAction", "Number/Dimension of actions per element", "1")

CALL prms%CreateRealOption("SR_reward_blendfac", "Exponential blending factor between [0,1] for reward per time step.\n"//&
                                                 "0.: Variable is not updated, 1.: always take new value.", "0.03")
CALL prms%CreateRealOption("SR_action_blendfac", "Exponential blending factor between [0,1] for actions per time step.\n"//&
                                                 "0.: Variable is not updated, 1.: always take new value.", "0.01")
CALL prms%CreateLogicalOption("SR_ClusteredDatabase", "SmartRedis database is clustered", ".FALSE.")

END SUBROUTINE DefineParametersSmartRedis


!==================================================================================================================================
!> Initialize SmartRedis client and allocate necessary tensors
!==================================================================================================================================
SUBROUTINE InitSmartRedis()
! MODULES
USE MOD_Globals
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars,       ONLY: nBCs,BC
USE MOD_ReadInTools,     ONLY: GETLOGICAL,GETREAL,GETINT,GETINTFROMSTR
USE MOD_TimeDisc_Vars,   ONLY: nRKStages
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: i
!==================================================================================================================================

doSmartRedis  = GETLOGICAL("doSmartRedis")
IF (doSmartRedis) THEN
  dbIsClustered = GETLOGICAL("SR_ClusteredDatabase")
  ! Currently only the MPI root communicates with the Database. Could be changing in the future.
  IF(MPIroot) SR_Error = Client%Initialize(dbIsClustered)
END IF

SR_Type = GETINTFROMSTR('SR_Type')

SELECT CASE (SR_Type)
CASE (PRM_SMARTREDIS_NONE)
  RETURN
CASE (PRM_SMARTREDIS_HIT)
  useInvariants = GETLOGICAL("SR_useInvariants")
  IF (useInvariants) doNormInvariants = GETLOGICAL("SR_doNormInvariants")
  SR_nVarAction  = GETINT("SR_nVarAction")
#if FV_ENABLED == 2
  IF (SR_nVarAction.NE.1) CALL ABORT(__STAMP__,
      'Only one action per element is supported for FV_ENABLED == 2')
#endif
#if EDDYVISCOSITY
  IF ( (SR_nVarAction.GT.2) .OR. (SR_nVarAction.LT.1) ) CALL ABORT(__STAMP__, &
      'Only one or two actions per element are supported for EDDYVISCOSITY')
#endif
CASE (PRM_SMARTREDIS_CHANNEL)
  CALL ABORT(__STAMP__, 'CHANNEL case not implemented yet')
CASE (PRM_SMARTREDIS_CYLINDER)
  SR_nVarAction  = GETINT("SR_nVarAction")
  IF (SR_nVarAction.NE.2) CALL ABORT(__STAMP__, &
      'Exactly two actions per element are supported for CYLINDER case')
  ! If we have cylinder, identify BC used for reward
  IF (SR_Type.EQ.PRM_SMARTREDIS_CYLINDER) THEN
    IF (COUNT(BC.EQ.31).NE.1) CALL ABORT(__STAMP__, &
        'Exactly one BC of type 31 (cylinder) must be defined for SmartRedis cylinder case')
    DO i=1,nBCs
      IF (BC(i).EQ.31) THEN
        SR_BC = i
        EXIT
      END IF
    END DO
    ! TODO: Compute automatically based on analyze_dt and relevant used-defined time scales
    SR_reward_blendfac = GETREAL("SR_reward_blendfac")
    SR_action_blendfac = GETREAL("SR_action_blendfac")
    IF (SR_reward_blendfac.LT.0. .OR. SR_reward_blendfac.GT.1.) CALL ABORT(__STAMP__, &
        'SR_reward_blendfac must be in [0,1]')
    IF (SR_action_blendfac.LT.0. .OR. SR_action_blendfac.GT.1.) CALL ABORT(__STAMP__, &
        'SR_action_blendfac must be in [0,1]')
    ! Nullify temporally filtered quantities
    SR_actions = 0.
    SR_BodyForce = 0.
  END IF
CASE DEFAULT
  CALL ABORT(__STAMP__, 'Unknown SmartRedis communication type')
END SELECT
END SUBROUTINE InitSmartRedis


!==================================================================================================================================
!> First gather arbitrary array RealArray across all ranks and communicate that to SmartRedis
!==================================================================================================================================
SUBROUTINE GatheredWriteSmartRedis(rank, nVal, RealArray, key, Shape_Out)
! MODULES
USE MOD_Globals
USE MOD_SmartRedis_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: rank                      !< Rank of array
INTEGER,INTENT(IN)             :: nVal(rank)                !< Local number of variables in each rank
REAL,INTENT(IN)                :: RealArray(PRODUCT(nVal))  !< Real array to write
CHARACTER(LEN=*),INTENT(IN)    :: Key                       !< array name to write to database
INTEGER,INTENT(IN),OPTIONAL    :: Shape_Out(rank)           !< shape of output array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: RealArray_Global(:)
INTEGER                        :: nValGather(rank)
#if USE_MPI
INTEGER                        :: i,nDOFLocal
INTEGER,DIMENSION(nProcessors) :: nDOFPerRank,offsetRank
#endif
!==================================================================================================================================
#if USE_MPI
nDOFLocal=PRODUCT(nVal)
CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerRank,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)

offsetRank=0
IF(MPIroot)THEN
  nValGather=nVal
  nValGather(rank)=SUM(nDOFPerRank)/PRODUCT(nVal(1:rank-1))
  DO i=2,nProcessors
    offsetRank(i)=offsetRank(i-1)+nDOFPerRank(i-1)
  END DO
  ALLOCATE(RealArray_Global(PRODUCT(nValGather)))
ELSE
  ALLOCATE(RealArray_Global(1))
ENDIF

CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                 RealArray_Global,nDOFPerRank,offsetRank,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#else
nValGather = nVal
ALLOCATE(RealArray_Global(PRODUCT(nValGather)))
RealArray_Global = RealArray
#endif /*USE_MPI*/

IF(MPIroot)THEN
  IF (PRESENT(Shape_Out)) THEN
    IF (PRODUCT(nValGather) .NE. PRODUCT(Shape_Out)) CALL ABORT(__STAMP__, &
                                                                'Wrong output dimension in GatheredWrite to SmartRedis!')
    SR_Error = Client%put_tensor(TRIM(Key), REAL(RealArray_Global,4), Shape_Out)
  ELSE
    SR_Error = Client%put_tensor(TRIM(Key), REAL(RealArray_Global,4), SHAPE(RealArray_Global))
  ENDIF
ENDIF

DEALLOCATE(RealArray_Global)

END SUBROUTINE GatheredWriteSmartRedis


!==================================================================================================================================
!> Get array from SmartRedis and Scatter it across all ranks
!==================================================================================================================================
SUBROUTINE GatheredReadSmartRedis(rank, nVal, RealArray, key)
! MODULES
USE MOD_Globals
USE MOD_SmartRedis_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: rank                      !< Rank of array
INTEGER,INTENT(IN)             :: nVal(rank)                !< Local number of variables in each rank
REAL,INTENT(INOUT)             :: RealArray(PRODUCT(nVal))  !< Real array to write
CHARACTER(LEN=*),INTENT(IN)    :: Key                       !< array name to write to database
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: RealArray_Global(:)
#if USE_MPI
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nProcessors) :: nDOFPerRank,offsetRank
#endif
INTEGER,PARAMETER              :: interval = 10      ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1) ! Infinite number of polling tries
LOGICAL                        :: found    = .FALSE.
!==================================================================================================================================
#if USE_MPI
nDOFLocal=PRODUCT(nVal)
CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerRank,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)

offsetRank=0
IF(MPIroot) THEN
  nValGather=nVal
  nValGather(rank)=SUM(nDOFPerRank)/PRODUCT(nVal(1:rank-1))
  DO i=2,nProcessors
    offsetRank(i)=offsetRank(i-1)+nDOFPerRank(i-1)
  END DO
  ALLOCATE(RealArray_Global(PRODUCT(nValGather)))
ELSE
  ALLOCATE(RealArray_Global(1))
ENDIF
#else
ALLOCATE(RealArray_Global(PRODUCT(nVal)))
#endif /*USE_MPI*/

IF(MPIroot) THEN
  SR_Error = Client%poll_tensor(TRIM(Key), interval, tries, found)
  IF(.NOT. found) CALL ABORT(__STAMP__, 'Failed to retrieve tensor with key '//TRIM(key))
  SR_Error = Client%unpack_tensor(TRIM(Key), RealArray_Global, SHAPE(RealArray_Global))
  SR_Error = Client%delete_tensor(TRIM(Key))
ENDIF

#if USE_MPI
CALL MPI_ScatterV(RealArray_Global,nDOFPerRank,offsetRank,MPI_DOUBLE_PRECISION,&
                  RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#else
RealArray = RealArray_Global
#endif
DEALLOCATE(RealArray_Global)

END SUBROUTINE GatheredReadSmartRedis

!==================================================================================================================================
!> Routine executed after each time step to pre/postprocess actions, states and rewards
!==================================================================================================================================
SUBROUTINE AnalyzeSmartRedis()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars,      ONLY: nBCs
USE MOD_CalcBodyForces, ONLY: CalcBodyForces
USE MOD_Exactfunc_Vars, ONLY: jetStrength
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(3,nBCs)    :: Fv,Fp,BodyForce ! viscous/pressure/resulting body force
!==================================================================================================================================
SELECT CASE (SR_Type)
CASE (PRM_SMARTREDIS_NONE)
  RETURN
CASE (PRM_SMARTREDIS_HIT)
  RETURN
CASE (PRM_SMARTREDIS_CHANNEL)
  RETURN
CASE (PRM_SMARTREDIS_CYLINDER)
  ! Update actions and bodyforces for reward using exponential filtering to avoid spikes.
  ! formula for exponential blending: y = y + alpha*(x-y) or y = alpha*x + (1-alpha)*y
  !
  ! 1.) Update Body Force (exponential blending)
  CALL CalcBodyForces(Bodyforce, Fp, Fv)
  SR_BodyForce = SR_BodyForce + SR_reward_blendfac*(BodyForce(:,SR_BC) - SR_BodyForce)
  ! 2.) Update Jet Strength (exponential blending)
  jetStrength = jetStrength + SR_action_blendfac*(SR_actions - jetStrength)
CASE DEFAULT
  CALL ABORT(__STAMP__, 'Unknown SmartRedis communication type')
END SELECT
END SUBROUTINE AnalyzeSmartRedis

!==================================================================================================================================
!> Call corresponding subroutine for the SmartRedis communication type
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis(U, firstTimeStep, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars, ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems)
LOGICAL,INTENT(IN)          :: FirstTimeStep
LOGICAL,INTENT(IN)          :: LastTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SELECT CASE (SR_Type)
CASE (PRM_SMARTREDIS_NONE)
  RETURN
CASE (PRM_SMARTREDIS_HIT)
  CALL ExchangeDataSmartRedis_HIT(U, FirstTimeStep, LastTimeStep)
CASE (PRM_SMARTREDIS_CHANNEL)
  CALL ExchangeDataSmartRedis_CHANNEL(U, FirstTimeStep, LastTimeStep)
CASE (PRM_SMARTREDIS_CYLINDER)
  CALL ExchangeDataSmartRedis_CYLINDER(U, FirstTimeStep, LastTimeStep)
CASE DEFAULT
  CALL ABORT(__STAMP__, 'Unknown SmartRedis communication type')
END SELECT
END SUBROUTINE ExchangeDataSmartRedis

!==================================================================================================================================
!> Data exchange with SmartRedis for the HIT case
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis_HIT(U, firstTimeStep, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_ChangeBasisByDim,   ONLY: ChangeBasisVolume
USE MOD_Interpolation_Vars, ONLY: Vdm_Leg
USE MOD_Mesh_Vars,          ONLY: nElems,nGlobalElems
USE MOD_Lifting_Vars,       ONLY: gradUx,gradUy,gradUz
#if USE_FFTW
USE MOD_FFT_Vars,           ONLY: kmax
USE MOD_Testcase_Vars,      ONLY: E_k
#endif
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,      ONLY: Cs
#endif
#if FV_ENABLED == 2
USE MOD_FV_Vars,            ONLY: FV_alpha,FV_alpha_min,FV_alpha_max
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems)
LOGICAL,INTENT(IN)          :: FirstTimeStep
LOGICAL,INTENT(IN)          :: LastTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: Key
REAL                           :: inv(5,0:PP_N,0:PP_N,0:PP_N,1:nElems)
REAL                           :: actions(SR_nVarAction,nElems)
REAL                           :: actions_modal(1,0:PP_N,0:PP_N,0:PP_N,nElems)
REAL                           :: Vdm(0:PP_N,0:PP_N)
INTEGER                        :: lastTimeStepInt(1),Dims(5),Dims_Out(5)
INTEGER                        :: i,iElem
#if EDDYVISCOSITY
INTEGER                        :: j,k
#endif
INTEGER,PARAMETER              :: interval = 10   ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1)   ! Infinite number of polling tries
!==================================================================================================================================
! Gather U across all MPI ranks and write to Redis Database
Key = TRIM(FlexiTag)//"state"
IF (useInvariants) THEN
  CALL ComputeInvariants(gradUx, gradUy, gradUz, inv, doNormInvariants)

  Dims = SHAPE(inv)
  Dims_Out(:) = Dims(:)
  Dims_Out(5) = nGlobalElems

  CALL GatheredWriteSmartRedis(5, Dims, inv, TRIM(Key), Shape_Out = Dims_Out)
ELSE
  Dims = SHAPE(U)
  Dims_Out(:) = Dims(:)
  Dims_Out(5) = nGlobalElems

  CALL GatheredWriteSmartRedis(5, Dims, U, TRIM(Key), Shape_Out = Dims_Out)
END IF

IF (MPIroot .AND. (.NOT. firstTimeStep)) THEN
#if USE_FFTW
  ! Put Energy Spectrum into DB for Reward
  Key = TRIM(FlexiTag)//"Ekin"
  SR_Error = Client%put_tensor(TRIM(Key),E_k,SHAPE(E_k))
#endif

  ! Indicate if FLEXI is about to finalize
  lastTimeStepInt = MERGE(-1,1,lastTimeStep)
  Key = TRIM(FlexiTag)//"step_type"
  SR_Error = Client%put_tensor(TRIM(Key),lastTimeStepInt,(/1/))
ENDIF

! Get Cs from Redis Database and scatter across all MPI ranks
! Only necessary if we want to compute further, i.e. if not lastTimeStep
IF (.NOT. lastTimeStep) THEN
  Key = TRIM(FlexiTag)//"actions"
  CALL GatheredReadSmartRedis(SIZE(SHAPE(actions)), SHAPE(actions), actions, TRIM(Key))

  ! Construct field data from obtained modes
  actions_modal = 0.
  actions_modal(1,0,0,0,:) = actions(1,:) ! constant
  IF (SR_nVarAction.EQ.2) THEN
    actions_modal(1,2,0,0,:) = actions(2,:)-0.25 ! quad. x-direction
    actions_modal(1,0,2,0,:) = actions(2,:)-0.25 ! quad. y-direction
    actions_modal(1,0,0,2,:) = actions(2,:)-0.25 ! quad. z-direction
  END IF

  ! Build the non-normalized VDM from modal to nodal
  ! Means all Legendre-Polynomials fulfill P(1) = 1 (makes coefficients more interpretable)
  DO i=0,PP_N
    Vdm(:,i) = Vdm_Leg(:,i)/SQRT(REAL(i)+0.5)
  END DO
  DO iElem=1,nElems
#if EDDYVISCOSITY
    CALL ChangeBasisVolume(PP_N,PP_N,Vdm,actions_modal(1,:,:,:,iElem),Cs(      1,:,:,:,iElem))
    ! Limit resulting Cs to specificed range
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Cs(1,i,j,k,iElem) = MAX(0., MIN(0.5, Cs(1,i,j,k,iElem)))
    END DO; END DO; END DO
#elif FV_ENABLED == 2
    ! Limit resulting alpha to specificed range
    FV_alpha(iElem) = MAX(FV_alpha_min, MIN(FV_alpha_max, actions_modal(1,0,0,0,iElem)))
#endif
  END DO
END IF !.NOT.lastTimeStep
END SUBROUTINE ExchangeDataSmartRedis_HIT

!==================================================================================================================================
!> Data exchange with SmartRedis for the CHANNEL case
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis_CHANNEL(U, firstTimeStep, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,          ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems)
LOGICAL,INTENT(IN)          :: FirstTimeStep
LOGICAL,INTENT(IN)          :: LastTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL ABORT(__STAMP__, 'CHANNEL case not implemented yet')
END SUBROUTINE ExchangeDataSmartRedis_CHANNEL

!==================================================================================================================================
!> Data exchange with SmartRedis for the CYLINDER case
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis_CYLINDER(U, firstTimeStep, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars,          ONLY: nElems
USE MOD_RecordPoints,       ONLY: EvalRecordPoints
USE MOD_RecordPoints_Vars,  ONLY: nRP,nGlobalRP
USE MOD_EOS,                ONLY: ConsToPrim
USE MOD_Exactfunc_Vars,     ONLY: jetStrength
USE MOD_Equation_Vars,      ONLY: RefStatePrim,IniRefState
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems)
LOGICAL,INTENT(IN)          :: FirstTimeStep
LOGICAL,INTENT(IN)          :: LastTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: Key
INTEGER,PARAMETER              :: interval = 10      ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1) ! Infinite number of polling tries
REAL                           :: U_RP(    PP_nVar    ,nRP)  ! cons. state at record points
REAL                           :: UPrim_RP(PP_nVarPrim,nRP)  ! prim. state at record points
LOGICAL                        :: found    = .FALSE.
INTEGER                        :: i,lastTimeStepInt
REAL                           :: actions(1),cd,cl
!==================================================================================================================================
! Gather U across all MPI ranks and write to Redis Database
Key = TRIM(FlexiTag)//"state"
CALL EvalRecordPoints(U_RP)
DO i=1,nRP
  CALL ConsToPrim(UPrim_RP(:,i),U_RP(:,i))
  UPrim_RP(PRES,i) = UPrim_RP(PRES,i) - RefStatePrim(PRES,IniRefState) ! Subtract mean pressure

END DO
CALL GatheredWriteSmartRedis(1, SHAPE(UPrim_RP(PRES,:)), UPrim_RP(PRES,:), TRIM(Key), Shape_Out = (/nGlobalRP/))

IF (MPIroot .AND. (.NOT. firstTimeStep)) THEN
  ! Compute lift and drag coefficients, i.e. (area already taken into account in computation of forces)
  ! cd = Drag/(0.5*rho_0*U_0^2)
  ! cl = Lift/(0.5*rho_0*U_0^2)
  ! since U_bulk = rho_0 = 1 in this case, it simplifies to
  cd = SR_BodyForce(1)/0.5
  cl = SR_BodyForce(2)/0.5
  ! Put lift and drag coeffcicients into DB for Reward
  Key = TRIM(FlexiTag)//"reward"
  SR_Error = Client%put_tensor(TRIM(Key),(/cd,cl/),(/2/))

  ! Indicate if FLEXI is about to finalize
  lastTimeStepInt = MERGE(-1,1,lastTimeStep)
  Key = TRIM(FlexiTag)//"step_type"
  SR_Error = Client%put_tensor(TRIM(Key),lastTimeStepInt,(/1/))
ENDIF

! Get actions from Redis Database and scatter across all MPI ranks
! Only necessary if we want to compute further, i.e. if not lastTimeStep
IF (.NOT. lastTimeStep) THEN
  Key = TRIM(FlexiTag)//"actions"
  IF(MPIroot) THEN
    SR_Error = Client%poll_tensor(TRIM(Key), interval, tries, found)
    IF(.NOT. found) CALL ABORT(__STAMP__, 'Failed to retrieve tensor with key '//TRIM(key))
    SR_Error = Client%unpack_tensor(TRIM(Key), actions, SHAPE(actions))
    SR_Error = Client%delete_tensor(TRIM(Key))
  ENDIF
  ! Ensure zero net massflow
  SR_actions(1) =     actions(1)
  SR_actions(2) = -1.*actions(1)
#if USE_MPI
  CALL MPI_BCAST(SR_actions,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
#endif
END IF !.NOT.lastTimeStep
END SUBROUTINE ExchangeDataSmartRedis_CYLINDER

!==================================================================================================================================
!> Computes the 5 invariants of the velocity gradient tensor, according to
!>   - 'A more general effective-viscosity hypothesis', Pope, J. Fluid Mech.,1975.
!> NOTE: The first of the 6 invariants (trace(S)) of 'matrix' is trivially zero for incompressible flows, thus resulting in only
!>       5 non-trivial invariants.
!==================================================================================================================================
PPURE SUBROUTINE ComputeInvariants(gradUx, gradUy, gradUz, Invariants, doNormalize)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars ,ONLY: nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: gradUx(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(IN)             :: gradUy(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL,INTENT(IN)             :: gradUz(PP_nVarLifting,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
LOGICAL,INTENT(IN)          :: doNormalize
REAL,INTENT(INOUT)          :: Invariants(5,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: i,j,k,iElem
REAL                   :: S(3,3),S2(3,3)
REAL                   :: W(3,3),W2(3,3)
REAL                   :: mat(3,3)
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_NZ;DO j=0,PP_N; DO i=0,PP_N
    ! Velocity gradient tensor \nabla u
    mat(:,1) = gradUx(VELV,i,j,k,iElem)
    mat(:,2) = gradUy(VELV,i,j,k,iElem)
    mat(:,3) = gradUz(VELV,i,j,k,iElem)

    S = 0.5*(mat+TRANSPOSE(mat)) ! Symmetric part
    W = 0.5*(mat-TRANSPOSE(mat)) ! Anti-Symmetric part

    ! Normalize and de-dimensionalize invariants with the Frobenius norm of |S|=\sqrt{2*S_ij*S_ij}.
    IF (doNormalize) THEN
      ! This corresponds to making the velocity gradient tensor dimensionless with the characteristic velocity and length scale.
      ! This again corresponds to using the resolved dissipation rate and the viscosity.
      W = 1./SQRT(2*SUM(S(:,:)**2)) * W
      S = 1./SQRT(2*SUM(S(:,:)**2)) * S
    END IF

    ! Trace(S^2)
    S2 = MATMUL(S,S)
    Invariants(1,i,j,k,iElem) =  S2(1,1) +  S2(2,2) +  S2(3,3)

    ! Trace(W^2)
    W2 = MATMUL(W,W)
    Invariants(2,i,j,k,iElem) =  W2(1,1) +  W2(2,2) +  W2(3,3)

    ! Trace(S^3)
    mat = MATMUL(S2,S)
    Invariants(3,i,j,k,iElem) = mat(1,1) + mat(2,2) + mat(3,3)

    ! Trace(W^2*S)
    mat = MATMUL(W2,S)
    Invariants(4,i,j,k,iElem) = mat(1,1) + mat(2,2) + mat(3,3)

    ! Trace(W^2*S^2)
    mat = MATMUL(W2,S2)
    Invariants(5,i,j,k,iElem) = mat(1,1) + mat(2,2) + mat(3,3)

  END DO; END DO; END DO
END DO
END SUBROUTINE ComputeInvariants

!==================================================================================================================================
!> Finalize SmartRedis Client
!==================================================================================================================================
SUBROUTINE FinalizeSmartRedis()
! MODULES
USE MOD_Globals
USE MOD_SmartRedis_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
IF(MPIroot.AND.doSmartRedis) SR_Error = Client%destructor()
END SUBROUTINE FinalizeSmartRedis
#endif /*USE_SMARTREDIS*/

END MODULE MOD_SmartRedis
