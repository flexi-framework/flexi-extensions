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

REAL,ALLOCATABLE,TARGET            :: BodyForces(:)
REAL,ALLOCATABLE,TARGET            :: BodyForcesFine(:)
REAL,ALLOCATABLE,TARGET            :: BodyForcesFineSum(:)
REAL,ALLOCATABLE,TARGET            :: BodyForcesFineSqSum(:)

REAL,ALLOCATABLE,TARGET            :: BodyForcesCoarse(:)
REAL,ALLOCATABLE,TARGET            :: BodyForcesCoarseSum(:)
REAL,ALLOCATABLE,TARGET            :: BodyForcesCoarseSqSum(:)

REAL,ALLOCATABLE,TARGET            :: DBodyForcesSqSum(:)
REAL,ALLOCATABLE,TARGET            :: SigmaSqBodyForces(:)
REAL,ALLOCATABLE,TARGET            :: BiasBodyForces(:)

REAL,ALLOCATABLE,TARGET            :: MeanBodyForces(:)
REAL,ALLOCATABLE,TARGET            :: StdDevBodyForces(:)
INTEGER                            :: nValBodyForces(1)

CHARACTER(LEN=255)                 :: StateFileFine
CHARACTER(LEN=255)                 :: BodyFocesFileFine
CHARACTER(LEN=255)                 :: StateFileCoarse
CHARACTER(LEN=255)                 :: BodyFocesFileCoarse
CHARACTER(LEN=255)                 :: FileNameSums
CHARACTER(LEN=255)                 :: FileNameSumsBF
CHARACTER(LEN=255)                 :: FilenameMean
CHARACTER(LEN=255)                 :: FilenameVariance
CHARACTER(LEN=255)                 :: suffix
CHARACTER(LEN=255)                 :: FileType
CHARACTER(LEN=255)                 :: DataSetName
CHARACTER(LEN=255)                 :: DataSetNameBodyForces
CHARACTER(LEN=255)                 :: MeshFile_Sums
REAL                               :: Time_Sums


INTEGER                            :: nVarTotal=10
CHARACTER(LEN=255)                 :: VarNames(10)

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
