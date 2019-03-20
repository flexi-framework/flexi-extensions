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
CHARACTER(LEN=255)                 :: FileNameTS 
CHARACTER(LEN=255)                 :: FileNameFFT 

INTEGER                            :: RP_specified

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

LOGICAL                            :: doSPL
END MODULE MOD_Nisp_RP_Vars
