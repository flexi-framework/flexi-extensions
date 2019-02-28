#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_EstimateSigma_RP_Output
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteSumsToHDF5
  MODULE PROCEDURE WriteSumsToHDF5
END INTERFACE

INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

INTERFACE WriteSigmaSq
  MODULE PROCEDURE WriteSigmaSq
END INTERFACE


PUBLIC:: WriteSumsToHDF5
PUBLIC:: WriteMeanAndVarianceToHDF5
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write the sum of all deltaU and deltaU^2 RP data to HDF5 format for a specififc level
!===================================================================================================================================
SUBROUTINE WriteSumsToHDF5()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars   ,ONLY:  nPoints
USE MOD_EstimateSigma_RP_Vars
USE MOD_spec_Vars             ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu        ,ONLY: nVarVisu,VarNameVisu
USE MOD_OutputRPVisu_HDF5
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i
REAL,ALLOCATABLE                 :: tmp(:,:,:)
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:)
!===================================================================================================================================
! Output spectra
ALLOCATE (tmp(5*nVarVisu,nPoints,nSamples_spec))
ALLOCATE (VarNames_tmp(5*nVarVisu))
DO i=1,nVarVisu
  VarNames_tmp(i)=TRIM(VarNameVisu(i))//'_F'
  VarNames_tmp(i+1*nVarVisu)=TRIM(VarNameVisu(i))//'_C'
  VarNames_tmp(i+2*nVarVisu)=TRIM(VarNameVisu(i))//'_FSq'
  VarNames_tmp(i+3*nVarVisu)=TRIM(VarNameVisu(i))//'_CSq'
  VarNames_tmp(i+4*nVarVisu)=TRIM(VarNameVisu(i))//'_DSq'
END DO

tmp(1            :  nVarVisu,:,:) = UFineSum
tmp(  nVarVisu+1:2*nVarVisu,:,:) = UCoarseSum
tmp(2*nVarVisu+1:3*nVarVisu,:,:) = UFineSqSum
tmp(3*nVarVisu+1:4*nVarVisu,:,:) = UCoarseSqSum
tmp(4*nVarVisu+1:           ,:,:) = DUSqSum

WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING SUMS OF SPECTRA TO ', FileNameSumsOut
  CALL WriteDataToHDF5(nSamples_spec,nPoints,5*nVarVisu,VarNames_tmp,RPData_freq,tmp,FileNameSumsOut)
WRITE(UNIT_StdOut,'(132("-"))')
CALL WriteAttribute(File_ID,'SigmaSq',1,RealScalar=SigmaSq)
CALL WriteAttribute(File_ID,'SigmaSqFine',1,RealScalar=SigmaSqFine)
CALL WriteAttribute(File_ID,'Bias',1,RealScalar=Bias)


SDEALLOCATE (tmp)
SDEALLOCATE (VarNames_tmp)
END SUBROUTINE WriteSumsToHDF5


!===================================================================================================================================
!> Subroutine to write mean and variance of the RP data spectrum to HDF5 format for a specififc level
!===================================================================================================================================
SUBROUTINE WriteMeanAndVarianceToHDF5()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars   ,ONLY:  nPoints
USE MOD_EstimateSigma_RP_Vars
USE MOD_spec_Vars             ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu        ,ONLY: nVarVisu,VarNameVisu
USE MOD_OutputRPVisu_HDF5
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:)
REAL,ALLOCATABLE                 :: tmp(:,:,:)
!===================================================================================================================================
! Output spectra
!CoordNames(1)='Frequency'
ALLOCATE (VarNames_tmp(nVarVisu))
ALLOCATE (tmp(nVarVisu,nPoints,nSamples_spec))
DO i=1,nVarVisu
  VarNames_tmp(i)=TRIM(VarNameVisu(i))
END DO
tmp(:,:,:)=snSamples*(UFineSum-UCoarseSum)

WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN OF SPECTRA TO ', FileNameMean
  CALL WriteDataToHDF5(nSamples_spec,nPoints,nVarVisu,VarNames_tmp,RPData_freq,tmp,FileNameMean)
WRITE(UNIT_StdOut,'(132("-"))')

tmp(:,:,:)=SigmaSqSpec

WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING VARIANCE OF SPECTRA TO ', FileNameVariance
  CALL WriteDataToHDF5(nSamples_spec,nPoints,nVarVisu,VarNames_tmp,RPData_freq,tmp,FileNameVariance)
WRITE(UNIT_StdOut,'(132("-"))')

SDEALLOCATE (VarNames_tmp)
SDEALLOCATE (tmp)

END SUBROUTINE WriteMeanAndVarianceToHDF5

END MODULE MOD_EstimateSigma_RP_Output
