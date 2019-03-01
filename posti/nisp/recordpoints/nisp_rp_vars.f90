MODULE MOD_Nisp_RP_Vars
!===================================================================================================================================
! Contains global variables provided by the post processing routines
!===================================================================================================================================
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------

REAL,ALLOCATABLE                   :: UMeanTimeseries(:,:,:), UVarTimeseries(:,:,:), UModeTimeseries(:,:,:) 
REAL,ALLOCATABLE                   :: UFFT(:,:,:,:),UMeanFFT(:,:,:), UVarFFT(:,:,:), UModeFFT(:,:,:), time(:)
REAL,ALLOCATABLE                   :: UTimeseries(:,:,:,:)
CHARACTER(LEN=255)                 ::  FileNameMeanTS = 'mean_timeseries.h5' ,FileNameVarianceTS = 'variance_timeseries.h5'
CHARACTER(LEN=255)                 ::  FileNameMeanFFT = 'mean_spec.h5' ,FileNameVarianceFFT = 'variance_spec.h5' 

INTEGER                             :: RP_specified

INTEGER                            :: M 
INTEGER                            :: nStochSamples 
INTEGER                            :: nStochVars 
CHARACTER(LEN=255),ALLOCATABLE     :: StochVarNames(:), Distributions(:)
REAL, ALLOCATABLE                  :: StochPoints(:,:)
REAL, ALLOCATABLE                  :: StochWeights(:)
REAL, ALLOCATABLE                  :: DistributionProps(:,:)

INTEGER                            :: nRPFiles 
CHARACTER(LEN=255),ALLOCATABLE     :: RPFiles(:)

REAL, ALLOCATABLE                  :: HermiteCoeffs(:,:)
INTEGER                            :: nStochCoeffs 
INTEGER,ALLOCATABLE                :: MultiIndex(:,:)
END MODULE MOD_Nisp_RP_Vars
