MODULE MOD_MC_Vars
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
INTEGER                            :: nStateFiles
INTEGER                            :: nSamples
INTEGER                            :: iFileSamples

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
INTEGER                            :: addVars
INTEGER                            :: nElemsNew
CHARACTER(LEN=255)                 :: Mesh
CHARACTER(LEN=255)                 :: NodeType
REAL                               :: Time_State
REAL                               :: RefState(5)
REAL, ALLOCATABLE                  :: HermiteCoeff(:,:)
REAL,ALLOCATABLE                   :: UPrim(:,:,:,:,:)
REAL,ALLOCATABLE                   :: UAddVars(:,:,:,:,:)
REAL,ALLOCATABLE                   :: UMean(:,:,:,:,:), UVar(:,:,:,:,:), UMeanSquare(:,:,:,:,:) ,U(:,:,:,:,:)
REAL,ALLOCATABLE                   :: USampleAdditionalVars(:,:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: ElemData(:,:,:)
REAL,ALLOCATABLE,TARGET            :: ElemData_loc(:,:)
INTEGER                            :: nVal_ElemData(15)
CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesAdditional(:)

CHARACTER(LEN=255)                 :: StateFile, StateFileName, ProjectName
CHARACTER(LEN=255)                 :: FileNameMean = 'Mean.h5' ,FileNameVariance = 'Variance.h5'
END MODULE MOD_MC_Vars
