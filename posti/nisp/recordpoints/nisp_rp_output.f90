!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
#include "flexi.h"

MODULE MOD_Nisp_RP_Output
IMPLICIT NONE
PRIVATE

INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

INTERFACE WriteRPToHDF5
  MODULE PROCEDURE WriteRPToHDF5
END INTERFACE

PUBLIC:: WriteMeanAndVarianceToHDF5, WriteRPToHDF5
CONTAINS

!===================================================================================================================================
!> Subroutine to write mean and variance of the RP data spectrum to HDF5 format for a specififc level
!===================================================================================================================================
SUBROUTINE WriteMeanAndVarianceToHDF5()
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars     
USE MOD_RPSetVisuVisu_Vars  ,ONLY:  nPoints
USE MOD_RPSetVisuVisu_Vars                        
USE MOD_OutputRPVisu_Vars
USE MOD_RPSetVisu                        
USE MOD_OutputRPVisu 
USE MOD_Nisp_RP_Vars        ,ONLY: time, FileNameMeanTS,FileNameVarianceTS,FileNameMeanFFT,FileNameVarianceFFT, UMeanTimeseries,UMeanFFT, UVarTimeseries, UVarFFT
USE MOD_spec_Vars           ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu      ,ONLY: nVarVisu,VarNameVisu
USE MOD_OutputRPVisu_HDF5   ,ONLY: WriteDataToHDF5
USE MOD_HDF5_Output         ,ONLY: WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_Output_Vars         ,ONLY: ProgramName,FileVersion,ProjectName
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
  VarNames_tmp(i)=TRIM(VarNameVisu(i))//'_Mean'
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
! MEAN: 
WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN OF TIMESERIES TO ', FileNameMeanTS
  CALL WriteDataToHDF5(nSamples_out,nPoints,nVarVisu,VarNames_tmp,time,UMeanTimeseries,FileNameMeanTS)
WRITE(UNIT_StdOut,'(132("-"))')


WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN OF SPECTRA TO ', FileNameMeanFFT
  CALL WriteDataToHDF5(nSamples_spec,nPoints,nVarVisu,VarNames_tmp,RPData_freq,UMeanFFT,FileNameMeanFFT)
WRITE(UNIT_StdOut,'(132("-"))')



DO i=1,nVarVisu
  VarNames_tmp(i)=TRIM(VarNameVisu(i))//'_Variance'
END DO

WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING VARIANCE OF TIMESERIES TO ', FileNameVarianceTS
  CALL WriteDataToHDF5(nSamples_out,nPoints,nVarVisu,VarNames_tmp,time,UVarTimeseries,FileNameVarianceTS)
WRITE(UNIT_StdOut,'(132("-"))')

WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING VARIANCE OF SPECTRA TO ', FileNameVarianceFFT
  CALL WriteDataToHDF5(nSamples_spec,nPoints,nVarVisu,VarNames_tmp,RPData_freq,UVarFFT,FileNameVarianceFFT)
WRITE(UNIT_StdOut,'(132("-"))')

SDEALLOCATE (VarNames_tmp)

END SUBROUTINE WriteMeanAndVarianceToHDF5

!===================================================================================================================================
!> Subroutine to write mean and variance of the RP data spectrum to HDF5 format for a specififc level
!===================================================================================================================================
SUBROUTINE WriteRPToHDF5(UFFT)
! MODULES
USE MOD_Globals
USE MOD_RPData_Vars     

USE MOD_RPSetVisuVisu_Vars                        
USE MOD_OutputRPVisu_Vars
USE MOD_RPSetVisu                        
USE MOD_OutputRPVisu 
USE MOD_Nisp_RP_Vars         ,ONLY: time, FileNameMeanTS,FileNameVarianceTS,FileNameMeanFFT,FileNameVarianceFFT, UMeanTimeseries,UMeanFFT, UVarTimeseries, UVarFFT
USE MOD_spec_Vars             ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu        ,ONLY: nVarVisu,VarNameVisu
USE MOD_OutputRPVisu_HDF5     ,ONLY: WriteDataToHDF5
USE MOD_HDF5_Output         ,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
REAL, INTENT(IN)                 :: UFFT(1:nVarVisu,nRP_global,nSamples_spec)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:), Filename
REAL,ALLOCATABLE                 :: tmp(:,:,:)
!===================================================================================================================================
! Output spectra
!CoordNames(1)='Frequency'
ALLOCATE (VarNames_tmp(nVarVisu))
ALLOCATE (tmp(nVarVisu,nPoints,nSamples_spec))
DO i=1,nVarVisu
  VarNames_tmp(i)=TRIM(VarNameVisu(i))//'_FFT'
END DO
!-----------------------------------------------------------------------------------------------------------------------------------
Filename= 'UFFT'
WRITE(UNIT_StdOut,'(132("-"))')
  WRITE(UNIT_stdOut,'(A,A)')' WRITING FFT OF TIMESERIES TO ', Filename
  CALL WriteDataToHDF5(nSamples_out,nPoints,nVarVisu,VarNames_tmp,time,UMeanTimeseries,FileNameMeanTS)
WRITE(UNIT_StdOut,'(132("-"))')


END SUBROUTINE WriteRPToHDF5

END MODULE MOD_Nisp_RP_Output
