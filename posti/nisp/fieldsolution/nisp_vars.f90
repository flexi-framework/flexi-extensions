MODULE MOD_Nisp_Vars
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
INTEGER                            :: M
INTEGER                            :: nStochSamples
INTEGER                            :: nStochVars
CHARACTER(LEN=255),ALLOCATABLE     :: StochVarNames(:), Distributions(:)
REAL, ALLOCATABLE                  :: StochPoints(:,:)
REAL, ALLOCATABLE                  :: StochWeights(:)
REAL, ALLOCATABLE                  :: DistributionProps(:,:)

REAL, ALLOCATABLE                  :: HermiteCoeffs(:,:)
INTEGER                            :: nStochCoeffs
INTEGER,ALLOCATABLE                :: MultiIndex(:,:)

INTEGER                            :: NNew
INTEGER                            :: nVar     ! GP-wise vars in HDF5
INTEGER                            :: nElemsNew
CHARACTER(LEN=255)                 :: Mesh
CHARACTER(LEN=255)                 :: NodeType
CHARACTER(LEN=255), ALLOCATABLE    :: VarNamesField(:), VarNames(:)
REAL                               :: Time_State
REAL, ALLOCATABLE                  :: HermiteCoeff(:,:)
REAL,ALLOCATABLE                   :: UPrim(:,:,:,:,:)
REAL,ALLOCATABLE                   :: UMean(:,:,:,:,:), UVar(:,:,:,:,:), UMode(:,:,:,:,:) ,U(:,:,:,:,:)
REAL,ALLOCATABLE                   :: USample(:,:,:,:,:,:)

CHARACTER(LEN=255)                 :: StateFile, StateFileName, ProjectName
CHARACTER(LEN=255)                 :: FileNameMean = 'Mean.h5' ,FileNameVariance = 'Variance.h5'
END MODULE MOD_Nisp_Vars
