#include "flexi.h"

MODULE MOD_SmartRedis

#if USE_SMARTREDIS
! MODULES
IMPLICIT NONE
PRIVATE

PUBLIC :: DefineParametersSmartRedis
PUBLIC :: InitSmartRedis
PUBLIC :: FinalizeSmartRedis
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

SUBROUTINE InitSmartRedis()
! MODULES
USE MOD_PreProc
USE MOD_SmartRedis_Vars
USE MOD_ReadInTools         ,ONLY:GETLOGICAL
USE MOD_Mesh_Vars           ,ONLY:nElems
IMPLICIT NONE

dbIsClustered = GETLOGICAL("ClusteredDatabase")
ALLOCATE(u_tensor(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(result_tensor(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

call client%initialize(dbIsClustered)

END SUBROUTINE InitSmartRedis

!==================================================================================================================================

SUBROUTINE FinalizeSmartRedis()
USE MOD_SmartRedis_Vars
IMPLICIT NONE

SDEALLOCATE(u_tensor)
SDEALLOCATE(result_tensor)

END SUBROUTINE FinalizeSmartRedis
#endif

END MODULE MOD_SmartRedis
