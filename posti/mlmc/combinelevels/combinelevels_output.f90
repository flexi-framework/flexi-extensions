#include "flexi.h"

MODULE MOD_Combinelevels_Output
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
USE MOD_IO_HDF5
USE MOD_HDF5_Output,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_CombineLevels_Vars
USE MOD_SwapMesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!CHARACTER(LEN=255)            :: VarNames(5),VarNamesAddField(5)
CHARACTER(LEN=255)            :: FileNameMean 
CHARACTER(LEN=255)            :: FileNameVariance 
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
! MEAN: 
FileNameMean      = 'mean_state.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND VARIANCE TO HDF5 FILE "',TRIM(FileNameMean),'" ... \n'

CALL GenerateFileSkeleton(TRIM(FileNameMean),&
                          'State',&
                          nVar_State,&
                          NNew,&
                          VarNamesMean,&
                          MeshFileNew,&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)
CALL OpenDataFile(TRIM(FileNameMean),.FALSE.,.TRUE.,readOnly=.FALSE.)
CALL WriteArray('DG_Solution',5,(/nVar_State,NNew+1,NNew+1,NNew+1,nElemsNew/),(/nVar_State,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    (/0,0,0,0,0/),.FALSE.,RealArray=UMean)
CALL CloseDataFile()


!-----------------------------------------------------------------------------------------------------------------------------------
! VARIANCE: 
FileNameVariance  = 'variance_state.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND VARIANCE TO HDF5 FILE "',TRIM(FileNameVariance),'" ...\n'

CALL GenerateFileSkeleton(TRIM(FileNameVariance),&
                          'State',&
                          nVar_State,&
                          NNew,&
                          VarNamesVariance,&
                          MeshFileNew,&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)
CALL OpenDataFile(TRIM(FileNameVariance),.FALSE.,.TRUE.,readOnly=.FALSE.)
CALL WriteArray('DG_Solution',5,(/nVar_State,NNew+1,NNew+1,NNew+1,nElemsNew/),(/nVar_State,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    !(/0,0,0,0,0/),.FALSE.,RealArray=SQRT(ABS(UVariance)))
    (/0,0,0,0,0/),.FALSE.,RealArray=UVariance)
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5


END MODULE MOD_Combinelevels_Output
