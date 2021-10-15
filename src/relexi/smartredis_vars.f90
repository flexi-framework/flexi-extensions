#include "flexi.h"

MODULE MOD_SmartRedis_Vars
#if USE_SMARTREDIS
! MODULES
USE iso_c_binding
USE MOD_PreProc
USE smartredis_client, ONLY: client_type
USE MOD_Mesh_Vars    , ONLY: nElems
IMPLICIT NONE
PUBLIC
SAVE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
type(client_type)                                       :: client
real(kind=c_double), dimension(:,:,:,:,:), allocatable  :: u_tensor
real, dimension(:,:,:,:,:), allocatable                 :: result_tensor
character(len=255)                                      :: u_tensor_key
logical                                                 :: dbIsClustered
#endif

END MODULE MOD_SmartRedis_Vars
