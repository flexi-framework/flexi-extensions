#include "flexi.h"

MODULE MOD_SmartRedis

#if USE_SMARTREDIS
! MODULES
IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersSmartRedis
  MODULE PROCEDURE DefineParametersSmartRedis
END INTERFACE

INTERFACE InitSmartRedis
  MODULE PROCEDURE InitSmartRedis
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
PUBLIC :: GatheredWriteSmartRedis
PUBLIC :: ExchangeDataSmartRedis
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters for SmartRedis
!==================================================================================================================================
SUBROUTINE DefineParametersSmartRedis()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("Relexi")
CALL prms%CreateLogicalOption("ClusteredDatabase", "SmartRedis database is clustered", value=".false.")

END SUBROUTINE DefineParametersSmartRedis


!==================================================================================================================================
!> Initialize SmartRedis client and allocate necessary tensors
!==================================================================================================================================
SUBROUTINE InitSmartRedis()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_ReadInTools         ,ONLY:GETLOGICAL
USE MOD_Mesh_Vars           ,ONLY:nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

dbIsClustered = GETLOGICAL("ClusteredDatabase")

! Currently only the MPI root communicates with the Database. Could be changing in the future.
IF(MPIroot) CALL Client%Initialize(dbIsClustered)

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
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nProcessors) :: nDOFPerRank,offsetRank
!==================================================================================================================================

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

IF(MPIroot)THEN
  IF (PRESENT(Shape_Out)) THEN
    IF (PRODUCT(nValGather) .NE. PRODUCT(Shape_Out)) CALL ABORT(__STAMP__, &
                                                                'Wrong output dimension in GatheredWrite to SmartRedis!')
    CALL Client%put_tensor(TRIM(Key), REAL(RealArray_Global,4), Shape_Out)
  ELSE
    CALL Client%put_tensor(TRIM(Key), REAL(RealArray_Global,4), SHAPE(RealArray_Global))
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
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nProcessors) :: nDOFPerRank,offsetRank
INTEGER,PARAMETER              :: interval = 10   ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1)   ! Infinite number of polling tries
LOGICAL                        :: found = .FALSE.
!==================================================================================================================================

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

IF(MPIroot) THEN
  found = Client%poll_tensor(TRIM(Key), interval, tries)
  !IF(.NOT. found) CALL ABORT(__STAMP__, 'Failed to retrieve tensor with key '//TRIM(key))
  CALL Client%unpack_tensor(TRIM(Key), RealArray_Global, SHAPE(RealArray_Global))
  CALL Client%delete_tensor(TRIM(Key))
ENDIF

CALL MPI_ScatterV(RealArray_Global,nDOFPerRank,offsetRank,MPI_DOUBLE_PRECISION,&
                  RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)

DEALLOCATE(RealArray_Global)

END SUBROUTINE GatheredReadSmartRedis


!==================================================================================================================================
!> 
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis(U, firstTimeStep, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars       ,ONLY: nElems,nGlobalElems
#if USE_FFTW
USE MOD_FFT_Vars,           ONLY: kmax
USE MOD_Testcase_Vars,      ONLY: E_k
#endif
#if EDDYVISCOSITY
USE MOD_EddyVisc_Vars,      ONLY: Cs
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
INTEGER                        :: lastTimeStepInt(1),Dims(5),Dims_Out(5)
INTEGER,PARAMETER              :: interval = 10   ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1)   ! Infinite number of polling tries
!==================================================================================================================================
! Gather U across all MPI ranks and write to Redis Database
Dims = SHAPE(U)
Dims_Out(:) = Dims(:)
Dims_Out(5) = nGlobalElems

Key = TRIM(FlexiTag)//"U"
CALL GatheredWriteSmartRedis(5, Dims, U, TRIM(Key), Shape_Out = Dims_Out)

IF (MPIroot .AND. (.NOT. firstTimeStep)) THEN
#if USE_FFTW
  ! Put Energy Spectrum into DB for Reward
  Key = TRIM(FlexiTag)//"Ekin"
  CALL Client%put_tensor(TRIM(Key),E_k(1:kmax),(/kmax/))
#endif

  ! Indicate if FLEXI is about to finalize
  lastTimeStepInt = MERGE(-1,1,lastTimeStep)
  Key = TRIM(FlexiTag)//"step_type"
  CALL Client%put_tensor(TRIM(Key),lastTimeStepInt,(/1/))
ENDIF

#if EDDYVISCOSITY
! Get Cs from Redis Database and scatter across all MPI ranks
! Only necessary if we want to compute further, i.e. if not lastTimeStep
IF (.NOT. lastTimeStep) THEN
  Key = TRIM(FlexiTag)//"Cs"
  CALL GatheredReadSmartRedis(SIZE(SHAPE(Cs)), SHAPE(Cs), Cs, TRIM(Key))
END IF
#endif /* EDDYVISCOSITY */

END SUBROUTINE ExchangeDataSmartRedis


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

IF(MPIroot) CALL Client%destructor()

END SUBROUTINE FinalizeSmartRedis
#endif

END MODULE MOD_SmartRedis
