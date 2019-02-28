MODULE MOD_EstimateSigma_RP_Vars
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

INTEGER                             :: nStart
INTEGER                             :: nEnd
CHARACTER(LEN=255)                  :: FileNameSums
CHARACTER(LEN=255)                  :: FilenameMean
CHARACTER(LEN=255)                  :: FilenameVariance

INTEGER                             :: iLevel
INTEGER                             :: nIter
INTEGER                             :: nIterIn
INTEGER                             :: VarAna

INTEGER                             :: RP_specified
REAL                                :: SigmaSq
REAL                                :: SigmaSqFine
REAL                                :: Bias
REAL,ALLOCATABLE,TARGET             :: SigmaSqSpec(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UFine(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UCoarse(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UFineSum(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UCoarseSum(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UFineSqSum(:,:,:)
REAL,ALLOCATABLE,TARGET             :: UCoarseSqSum(:,:,:)
REAL,ALLOCATABLE,TARGET             :: DUSqSum(:,:,:)

REAL                                :: snSamples
REAL                                :: snSamplesM1


END MODULE MOD_EstimateSigma_RP_Vars
