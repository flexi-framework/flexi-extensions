MODULE MOD_CombineLevels_RP_Vars
!===================================================================================================================================
! Contains global variables provided by the post processing routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! tate File Attribute --------------------------------------------------------------------------------------------------------------

CHARACTER(LEN=255)                  :: Filename
CHARACTER(LEN=255),ALLOCATABLE      :: VarNamesRP(:)

INTEGER                             :: nF
INTEGER                             :: nVal
INTEGER                             :: nIter
INTEGER                             :: nLevels

REAL,ALLOCATABLE,TARGET             :: RP_Freq(:)
REAL,ALLOCATABLE,TARGET             :: CoordinatesRP(:,:)
REAL,ALLOCATABLE,TARGET             :: Utmp(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UMean(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UVariance(:,:,:)

END MODULE MOD_CombineLevels_RP_Vars
