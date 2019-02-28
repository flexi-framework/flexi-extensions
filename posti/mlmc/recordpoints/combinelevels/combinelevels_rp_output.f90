#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_CombineLevels_RP_Output
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

PUBLIC:: WriteMeanAndVarianceToHDF5
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write mean and variance of the RP data spectrum to HDF5 format for a specififc level
!===================================================================================================================================
SUBROUTINE WriteMeanAndVarianceToHDF5()
! MODULES
USE MOD_Globals
USE MOD_RPSetVisuVisu_Vars   ,ONLY:  nPoints
USE MOD_CombineLevels_RP_Vars 
USE MOD_OutputRPVisu_HDF5
USE MOD_OutputRPVisu_Vars     ,ONLY: nCoords,CoordNames

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:)
!===================================================================================================================================
! Output spectra
!CoordNames(1)='Frequency'


nCoords=4
ALLOCATE (CoordNames(nCoords))
CoordNames(1)='Time'
CoordNames(2)='CoordinateX'
CoordNames(3)='CoordinateY'
CoordNames(4)='CoordinateZ'


ALLOCATE (VarNames_tmp(nVal))
DO i=1,nVal
  VarNames_tmp(i)=TRIM(VarNamesRP(i))//'_Mean'
END DO

FileNameMean='mean_spec.h5'
WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN OF SPECTRA TO ', FileNameMean
  CALL WriteDataToHDF5(nF,nPoints,nVal,VarNames_tmp,RP_freq,UMean,FileNameMean)
WRITE(UNIT_StdOut,'(132("-"))')

DO i=1,nVal
  VarNames_tmp(i)=TRIM(VarNamesRP(i))//'_Variance'
END DO

FileNameVariance='variance_spec.h5'
WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING VARIANCE OF SPECTRA TO ', FileNameVariance
  CALL WriteDataToHDF5(nF,nPoints,nVal,VarNames_tmp,RP_freq,UVariance,FileNameVariance)
WRITE(UNIT_StdOut,'(132("-"))')

SDEALLOCATE (VarNames_tmp)

END SUBROUTINE WriteMeanAndVarianceToHDF5


END MODULE MOD_CombineLevels_RP_Output
