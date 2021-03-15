MODULE MOD_MLMC_Vars
!===================================================================================================================================
! Contains global variables provided by the post processing routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
REAL,ALLOCATABLE,TARGET            :: UFine(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UFineSum(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UFineSqSum(:,:,:,:,:)

REAL,ALLOCATABLE,TARGET            :: UCoarse(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UCoarseSum(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UCoarseSqSum(:,:,:,:,:)

REAL,ALLOCATABLE,TARGET            :: DUSqSum(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: SigmaSqField(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: BiasField(:,:,:,:,:)

REAL,ALLOCATABLE,TARGET            :: Mean(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: StdDev(:,:,:,:,:)
INTEGER                            :: nVal(5)

CHARACTER(LEN=255)                 :: StateFileFine
CHARACTER(LEN=255)                 :: StateFileCoarse
CHARACTER(LEN=255)                 :: FileNameSums
CHARACTER(LEN=255)                 :: FilenameMean
CHARACTER(LEN=255)                 :: FilenameVariance
CHARACTER(LEN=255)                 :: suffix
CHARACTER(LEN=255)                 :: FileType
CHARACTER(LEN=255)                 :: DataSetName
CHARACTER(LEN=255)                 :: MeshFile_Sums
REAL                               :: Time_Sums

INTEGER                            :: nVarTotal=11
CHARACTER(LEN=255)                 :: VarNames(11)

INTEGER                            :: nStart
INTEGER                            :: nEnd
INTEGER                            :: nIter
INTEGER                            :: nIterIn
INTEGER                            :: nSamples_Sums
REAL                               :: snSamples
REAL                               :: snSamplesM1

INTEGER                            :: VarAna
INTEGER                            :: iSample
INTEGER                            :: iLevel
INTEGER                            :: nAna

REAL                               :: SigmaSq
REAL                               :: SigmaSqFine
REAL                               :: Bias

END MODULE MOD_MLMC_Vars
