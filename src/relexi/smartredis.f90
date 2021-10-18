#include "flexi.h"

MODULE MOD_SmartRedis

#if USE_SMARTREDIS
! MODULES
IMPLICIT NONE
PRIVATE

PUBLIC :: DefineParametersSmartRedis
PUBLIC :: InitSmartRedis
PUBLIC :: FinalizeSmartRedis
PUBLIC :: WriteToSmartRedis
!==================================================================================================================================

CONTAINS

!==================================================================================================================================

SUBROUTINE DefineParametersSmartRedis()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Relexi")
CALL prms%CreateLogicalOption("ClusteredDatabase", "SmartRedis database is clustered", value=".false.")

END SUBROUTINE DefineParametersSmartRedis
!==================================================================================================================================

SUBROUTINE WriteToSmartRedis(rank, nVal, RealArray)
! MODULES
USE MOD_Globals
USE MOD_SmartRedis_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)             :: rank                      !< Rank of array
INTEGER,INTENT(IN)             :: nVal(rank)                !< Local number of variables in each rank
REAL,INTENT(IN)                :: RealArray(PRODUCT(nVal))  !< Real array to write
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE               :: UReal(:)
INTEGER                        :: i,nValGather(rank),nDOFLocal
INTEGER,DIMENSION(nProcessors) :: nDOFPerNode,offsetNode
!==================================================================================================================================

nDOFLocal=PRODUCT(nVal)
CALL MPI_GATHER(nDOFLocal,1,MPI_INTEGER,nDOFPerNode,1,MPI_INTEGER,0,MPI_COMM_FLEXI,iError)

offsetNode=0
IF(myRank==0)THEN
    nValGather=nVal
    nValGather(rank)=SUM(nDOFPerNode)/PRODUCT(nVal(1:rank-1))
    DO i=2,nProcessors
        offsetNode(i)=offsetNode(i-1)+nDOFPerNode(i-1)
    END DO
    ALLOCATE(u_tensor(PRODUCT(nValGather)))
ELSE
    ALLOCATE(u_tensor(1))
ENDIF

CALL MPI_GATHERV(RealArray,nDOFLocal,MPI_DOUBLE_PRECISION,&
                 u_tensor,nDOFPerNode,offsetNode,MPI_DOUBLE_PRECISION,0,MPI_COMM_FLEXI,iError)

IF(myRank==0)THEN
    u_tensor_key = "u_tensor"

    call client%put_tensor(u_tensor_key, u_tensor, shape(u_tensor))
ENDIF

SDEALLOCATE(u_tensor)

end SUBROUTINE WriteToSmartRedis

!==================================================================================================================================

SUBROUTINE InitSmartRedis()
! MODULES
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_ReadInTools         ,ONLY:GETLOGICAL
USE MOD_Mesh_Vars           ,ONLY:nElems
IMPLICIT NONE

dbIsClustered = GETLOGICAL("ClusteredDatabase")
ALLOCATE(result_tensor(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

call client%initialize(dbIsClustered)

END SUBROUTINE InitSmartRedis

!==================================================================================================================================

SUBROUTINE FinalizeSmartRedis()
USE MOD_SmartRedis_Vars
IMPLICIT NONE

!SDEALLOCATE(u_tensor)
SDEALLOCATE(result_tensor)

END SUBROUTINE FinalizeSmartRedis
#endif

END MODULE MOD_SmartRedis
