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
INTEGER                            :: totSamples
INTEGER                            :: nStochVar
INTEGER                            :: NNew
INTEGER                            :: nVar     ! GP-wise vars in HDF5
INTEGER                            :: nElemsNew
INTEGER                            :: nPoly
INTEGER                            :: N0
CHARACTER(LEN=255)                 :: Mesh
CHARACTER(LEN=255)                 :: NodeType
CHARACTER(LEN=255), ALLOCATABLE    :: VarNamesField(:), VarNames(:)
REAL                               :: Time_State
REAL, ALLOCATABLE                  :: HermiteCoeff(:,:)
REAL, ALLOCATABLE                  :: weights(:)
REAL, ALLOCATABLE                  :: quadpoints(:,:), quadpointsTransformed(:,:),quadpointsInv(:,:)
INTEGER,ALLOCATABLE                :: MultiIndex(:,:)
REAL,ALLOCATABLE                   :: UPrim(:,:,:,:,:)
REAL,ALLOCATABLE                   :: UMean(:,:,:,:,:), UVar(:,:,:,:,:), UMode(:,:,:,:,:) ,U(:,:,:,:,:) 
INTEGER, ALLOCATABLE               :: distributions(:)

CHARACTER(LEN=255)                 :: StateFile, StateFileName
CHARACTER(LEN=255)                 :: sample_string, FileNameMean = 'Mean.h5' ,FileNameVariance = 'Variance.h5'  
CHARACTER(LEN=255)                 :: file_weights="parameter_weights.csv"
CHARACTER(LEN=255)                 :: file_quadpoints="parameter_quad.csv"
CHARACTER(LEN=255)                 :: file_quadpointsInv = "parameter_quadInv.csv"
CHARACTER(LEN=255)                 :: file_quadpointsTransformed="parameter_quadTransformed.csv"
END MODULE MOD_Nisp_Vars
