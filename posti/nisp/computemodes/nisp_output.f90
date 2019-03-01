#include "flexi.h"

MODULE MOD_Nisp_Output
!===================================================================================================================================
! Module for generic data output in HDF5 fromat
!===================================================================================================================================
! MODULES
USE HDF5
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

PUBLIC::WriteMeanAndVarianceToHDF5
!===================================================================================================================================
CONTAINS

SUBROUTINE WriteMeanAndVarianceToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_IO_HDF5
USE MOD_HDF5_Output         ,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_Nisp_Vars
USE MOD_Output_Vars,ONLY:ProgramName,FileVersion,ProjectName
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: Vars(2*nVar)
CHARACTER(LEN=255)            :: TimeString
Vars(1)='Density_Mean' 
Vars(2)='MomentumX_Mean' 
Vars(3)='MomentumY_Mean' 
Vars(4)='MomentumZ_Mean'
Vars(5)='EnergyStagnationDensity_Mean' 
Vars(6)='VelocityX_Mean' 
Vars(7)='VelocityY_Mean' 
Vars(8)='VelocityZ_Mean' 
Vars(9)='Pressure_Mean'
Vars(10)='Temperature_Mean'
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! MEAN: 
FileNameMean=TRIM(TIMESTAMP('Mean',Time_State))//'.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN TO HDF5 FILE "',TRIM(FileNameMean),'" ... \n'
CALL GenerateFileSkeleton(TRIM(FileNameMean),&
                          'State',&
                          2*nVar,&
                          NNew,&
                          Vars,&
                          'hopr_mesh.h5',&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)
                          
CALL OpenDataFile(TRIM(FileNameMean),.FALSE.,.TRUE.,readOnly=.FALSE.)
! WRITE DG SOLUTION ----------------------------------------------------------------------------------------------------------------
CALL WriteArray('DG_Solution',5,(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    (/0,0,0,0,0/),.FALSE.,RealArray=UMean(:,:,:,:,:))
CALL CloseDataFile()


! !-----------------------------------------------------------------------------------------------------------------------------------
! ! VARIANCE: 
Vars(1)='Density_Variance' 
Vars(2)='MomentumX_Variance' 
Vars(3)='MomentumY_Variance' 
Vars(4)='MomentumZ_Variance'
Vars(5)='EnergyStagnationDensity_Variance' 
Vars(6)='VelocityX_Variance' 
Vars(7)='VelocityY_Variance' 
Vars(8)='VelocityZ_Variance' 
Vars(9)='Pressure_Variance'
Vars(10)='Temperature_Variance'

FileNameVariance=TRIM(TIMESTAMP('Variance',Time_State))//'.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE VARIANCE TO HDF5 FILE "',TRIM(FileNameVariance),'" ... \n'

CALL GenerateFileSkeleton(TRIM(FileNameVariance),&
                          'State',&
                          2*nVar,&
                          NNew,&
                          Vars,&
                          'hopr_mesh.h5',&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)
! 
! ! Open HDF5 file for write
CALL OpenDataFile(TRIM(FileNameVariance),.FALSE.,.TRUE.,readOnly=.FALSE.)
! Write DG solution ----------------------------------------------------------------------------------------------------------------
CALL WriteArray('DG_Solution',5,(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    (/0,0,0,0,0/),.FALSE.,RealArray=UVar(:,:,:,:,:))
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5


END MODULE MOD_Nisp_Output
