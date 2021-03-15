#include "flexi.h"

MODULE MOD_MC_Output
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
INTERFACE WriteMeanAndVarianceMCToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceMCToHDF5
END INTERFACE

PUBLIC::WriteMeanAndVarianceMCToHDF5
!===================================================================================================================================
CONTAINS

SUBROUTINE WriteMeanAndVarianceMCToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Output         ,ONLY: WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
USE MOD_MC_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: VarNames(2*nVar+addVars)
CHARACTER(LEN=255)            :: FileName
VarNames( 1)='Density'
VarNames( 2)='MomentumX'
VarNames( 3)='MomentumY'
VarNames( 4)='MomentumZ'
VarNames( 5)='EnergyStagnationDensity'
VarNames( 6)='VelocityX'
VarNames( 7)='VelocityY'
VarNames( 8)='VelocityZ'
VarNames( 9)='Pressure'
VarNames(10)='Temperature'
VarNames(11)='cp'
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
FileName=TRIM(TIMESTAMP(TRIM(ProjectName)//'_'//TRIM('MC'),Time_State))//'.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND STDDEV TO HDF5 FILE "',TRIM(FileName),'" ... \n'
nGlobalElems=nElemsNew
CALL GenerateFileSkeleton(TRIM(FileName),'State',1,NNew,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                          mesh,Time_State,Time_State,withUserblock=.FALSE.,batchMode=.FALSE.,create=.TRUE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar+addVars,NNew,VarNames,&
                          mesh,Time_State,Time_State,create=.FALSE.,Dataset='Mean',batchMode=.FALSE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar+addVars,NNew,VarNames,&
                          mesh,Time_State,Time_State,create=.FALSE.,Dataset='StandardDeviation',batchMode=.FALSE.)

CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

CALL WriteArray('Mean',5,(/2*nVar+addVars,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar+addVars,NNew+1,NNew+1,NNew+1,nElemsNew/),(/0,0,0,0,0/),.FALSE.,RealArray=UMean)
CALL WriteArray('StandardDeviation',5,(/2*nVar+addVars,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar+addVars,NNew+1,NNew+1,NNew+1,nElemsNew/),(/0,0,0,0,0/),.FALSE.,RealArray=SQRT(UVar))
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceMCToHDF5


END MODULE MOD_MC_Output
