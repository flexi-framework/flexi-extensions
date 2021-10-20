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
ALLOCATE(result_tensor(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

CALL Client%Initialize(dbIsClustered)

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
INTEGER,DIMENSION(nProcessors) :: nDOFPerNode,offsetNode
!==================================================================================================================================

nDOFLocal=PRODUCT(nVal)
CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)

offsetNode=0
IF(MPIroot)THEN
  nValGather=nVal
  nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
  DO i=2,nProcessors
    offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
  END DO
  ALLOCATE(RealArray_Global(PRODUCT(nValGather)))
ELSE
  ALLOCATE(RealArray_Global(1))
ENDIF

CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                 RealArray_Global,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)

IF(MPIroot)THEN
  IF (PRESENT(Shape_Out)) THEN
    IF (PRODUCT(nValGather) .NE. PRODUCT(Shape_Out)) CALL ABORT(__STAMP__, &
                                                                'Wrong output dimension in GatheredWrite to SmartRedis!')
    CALL Client%put_tensor(TRIM(Key), RealArray_Global, Shape_Out)
  ELSE
    CALL Client%put_tensor(TRIM(Key), RealArray_Global, SHAPE(RealArray_Global))
  ENDIF
ENDIF

DEALLOCATE(RealArray_Global)

END SUBROUTINE GatheredWriteSmartRedis


!!==================================================================================================================================
!!> Get array from SmartRedis and Scatter it across all ranks
!!==================================================================================================================================
!SUBROUTINE GatheredReadSmartRedis(rank, nVal, RealArray, key, Shape_Out)
!! MODULES
!USE MOD_Globals
!USE MOD_SmartRedis_Vars
!IMPLICIT NONE
!!----------------------------------------------------------------------------------------------------------------------------------
!! INPUT/OUTPUT VARIABLES
!INTEGER,INTENT(IN)             :: rank                      !< Rank of array
!INTEGER,INTENT(IN)             :: nVal(rank)                !< Local number of variables in each rank
!REAL,INTENT(IN)                :: RealArray(PRODUCT(nVal))  !< Real array to write
!CHARACTER(LEN=*),INTENT(IN)    :: Key                       !< array name to write to database
!INTEGER,INTENT(IN),OPTIONAL    :: Shape_Out(rank)           !< shape of output array
!!----------------------------------------------------------------------------------------------------------------------------------
!! LOCAL VARIABLES
!REAL,ALLOCATABLE               :: RealArray_Global(:)
!INTEGER                        :: i,nValGather(rank),nDOFLocal
!INTEGER,DIMENSION(nProcessors) :: nDOFPerNode,offsetNode
!INTEGER,PARAMETER              :: interval = 10   ! polling interval in milliseconds
!INTEGER,PARAMETER              :: tries    = 1000 ! max. number of polling attempts before moving on
!!==================================================================================================================================
!
!nDOFLocal=PRODUCT(nVal)
!CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)
!
!offsetNode=0
!IF(MPIroot)THEN
!  nValGather=nVal
!  nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
!  DO i=2,nProcessors
!    offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
!  END DO
!  ALLOCATE(RealArray_Global(PRODUCT(nValGather)))
!ELSE
!  ALLOCATE(RealArray_Global(1))
!ENDIF
!
!!IF(MPIroot)THEN
!!  IF (PRESENT(Shape_Out)) THEN
!!    IF (PRODUCT(nValGather) .NE. PRODUCT(Shape_Out)) CALL ABORT(__STAMP__, &
!!                                                                'Wrong output dimension in GatheredWrite to SmartRedis!')
!!    CALL Client%put_tensor(TRIM(Key), RealArray_Global, Shape_Out)
!!  ELSE
!!    CALL Client%put_tensor(TRIM(Key), RealArray_Global, SHAPE(RealArray_Global))
!!  ENDIF
!!ENDIF
!IF(MPIroot)THEN
!  found = Client%poll_tensor(TRIM(Key), interval, tries)
!  CALL Client%unpack_tensor(TRIM(Key), RealArray, SHAPE(RealArray))
!  CALL Client%delete_tensor(TRIM(Key))
!ENDIF
!
!CALL MPI_ScatterV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
!                 RealArray_Global,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)
!
!DEALLOCATE(RealArray_Global)
!
!END SUBROUTINE GatheredReadSmartRedis


!==================================================================================================================================
!> 
!==================================================================================================================================
SUBROUTINE ExchangeDataSmartRedis(U, lastTimeStep)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_Mesh_Vars       ,ONLY: nElems,nGlobalElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)             :: U(1:3,0:PP_N,0:PP_N,0:PP_N,1:nElems)
LOGICAL,INTENT(IN)          :: LastTimeStep
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)             :: Key
LOGICAL                        :: found
INTEGER                        :: lastTimeStepInt(1),Dims(5),Dims_Out(5)
INTEGER,PARAMETER              :: interval = 10   ! polling interval in milliseconds
INTEGER,PARAMETER              :: tries    = HUGE(1) ! max. number of polling attempts before moving on
REAL                           :: Cs(nGlobalElems)
!==================================================================================================================================
Dims = SHAPE(U)
Dims_Out(:) = Dims(:)
Dims_Out(5) = nGlobalElems

Key = TRIM(FlexiTag)//"U"
CALL GatheredWriteSmartRedis(5, Dims, U, TRIM(Key), Shape_Out = Dims_Out)

IF (MPIroot) THEN
  lastTimeStepInt = MERGE(-1,1,lastTimeStep)
  Key = TRIM(FlexiTag)//"step_type"
  CALL Client%put_tensor(TRIM(Key),lastTimeStepInt,(/1/))

  Key = TRIM(FlexiTag)//"Cs"
  found = Client%poll_tensor(TRIM(Key), interval, tries)
  CALL Client%unpack_tensor(TRIM(Key), Cs, SHAPE(Cs))
  CALL Client%delete_tensor(TRIM(Key))

  ! TODO: Use ScatterV to scatter the received Cs to the corresponding MPI ranks.
  ! CALL GatheredReadSmartRedis(.....)

! TODO; Alibi to init Cs with smth on non-root ranks....
ELSE
  Cs = -1.
ENDIF

END SUBROUTINE ExchangeDataSmartRedis


!==================================================================================================================================
!> Finalize SmartRedis Client
!==================================================================================================================================
SUBROUTINE FinalizeSmartRedis()
! MODULES
USE MOD_SmartRedis_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

CALL Client%destructor()

END SUBROUTINE FinalizeSmartRedis
#endif

END MODULE MOD_SmartRedis
