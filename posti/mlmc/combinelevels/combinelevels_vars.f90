MODULE MOD_CombineLevels_Vars
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
CHARACTER(LEN=255)                 :: ProjectName(2)
CHARACTER(LEN=255)                 :: ProgramName(2) 
INTEGER                            :: N_HDF5(2)
INTEGER                            :: nVar_HDF5(2)     ! GP-wise vars in HDF5
INTEGER                            :: nGlobalElems_HDF5(2)
CHARACTER(LEN=255)                 :: NodeType_HDF5(2)
CHARACTER(LEN=255)                 :: MeshFile_HDF5(2)
REAL                               :: time_HDF5(2)
REAL,ALLOCATABLE                   :: Vdm(:,:)
!Attributes
CHARACTER(LEN=255)                 :: refinement,nAnaStr
CHARACTER(LEN=255)                 :: nLevelsStr
CHARACTER(LEN=255)                 :: VisuString
CHARACTER(LEN=255)                 :: FileType(2)
REAL                               :: FileVersion(2)
INTEGER                            :: nSamples(2)

INTEGER                            :: nAna,nAnaFinest,nAnaCoarsest,nElemsFinest,nElemsCoarsest




CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesMean(:)
CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesVariance(:)
CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesMean_Field(:)
CHARACTER(LEN=255),ALLOCATABLE     :: VarNamesVariance_Field(:)
INTEGER                            :: nLevel
INTEGER                            :: N_nLevel
REAL,ALLOCATABLE,TARGET            :: UMean(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UVariance(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UMeanField(:,:,:,:,:)
REAL,ALLOCATABLE,TARGET            :: UVarianceField(:,:,:,:,:)






END MODULE MOD_CombineLevels_Vars
