MODULE MOD_EstimateSigma_Vars
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

CHARACTER(LEN=255)                 :: StateFileFine
CHARACTER(LEN=255)                 :: StateFileCoarse
CHARACTER(LEN=255)                 :: FileNameSums
CHARACTER(LEN=255)                 :: FilenameMean
CHARACTER(LEN=255)                 :: FilenameVariance
CHARACTER(LEN=255)                 :: suffix
CHARACTER(LEN=255)                 :: FileType
CHARACTER(LEN=255)                 :: DataSetName

INTEGER                            :: nStart
INTEGER                            :: nEnd
INTEGER                            :: nIter
INTEGER                            :: nIterIn
REAL                               :: snSamples
REAL                               :: snSamplesM1

INTEGER                            :: VarAna
INTEGER                            :: nVarTotal
INTEGER                            :: iSample
INTEGER                            :: iLevel
INTEGER                            :: nAna

REAL                               :: SigmaSq
REAL                               :: SigmaSqFine
REAL                               :: Bias

INTEGER             :: FineNElemsOld                 !< Polynomial degree of old state
INTEGER             :: FineNState                    !< Polynomial degree of old state
!Vdm_GPNState_GPNNew needs to be DEALLOCATED
REAL,ALLOCATABLE    :: FineVdm_CLNInter_GPNNew(:,:)   !< Vandermonde from interpolation points on new mesh to solution points
REAL,ALLOCATABLE    :: FineVdm_GPNState_GPNNew(:,:)   !< Vandermonde from old solution to new solution (used for equal elements)
!ALL ALLOCATED IN INIT SWAMPMESH, ONLY DEPEND ON NEW VARS
REAL,ALLOCATABLE    :: FinexiInter(:,:,:,:,:)   !> Parametric coords of interpolation points of new solution in old mesh
INTEGER,ALLOCATABLE :: FineInterToElem(:,:,:,:) !> Mapping from new interpolation points  to elemID in old mesh
INTEGER,ALLOCATABLE :: FineequalElem(:)         !> Mapping between equal elements (iElemOld=equalElem(iElemNew))
LOGICAL,ALLOCATABLE :: FineIPDone(:,:,:,:)      !> Info if IP has been found


INTEGER             :: CoarseNElemsOld                 !< Polynomial degree of old state
INTEGER             :: CoarseNState                    !< Polynomial degree of old state
!Vdm_GPNState_GPNNew needs to be DEALLOCATED
REAL,ALLOCATABLE    :: CoarseVdm_CLNInter_GPNNew(:,:)   !< Vandermonde from interpolation points on new mesh to solution points
REAL,ALLOCATABLE    :: CoarseVdm_GPNState_GPNNew(:,:)   !< Vandermonde from old solution to new solution (used for equal elements)
!ALL ALLOCATED IN INIT SWAMPMESH, ONLY DEPEND ON NEW VARS
REAL,ALLOCATABLE    :: CoarsexiInter(:,:,:,:,:)   !> Parametric coords of interpolation points of new solution in old mesh
INTEGER,ALLOCATABLE :: CoarseInterToElem(:,:,:,:) !> Mapping from new interpolation points  to elemID in old mesh
INTEGER,ALLOCATABLE :: CoarseequalElem(:)         !> Mapping between equal elements (iElemOld=equalElem(iElemNew))
LOGICAL,ALLOCATABLE :: CoarseIPDone(:,:,:,:)      !> Info if IP has been found

END MODULE MOD_EstimateSigma_Vars
