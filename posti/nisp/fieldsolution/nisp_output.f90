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
USE MOD_IO_HDF5
USE MOD_HDF5_Output         ,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_Nisp_Vars
USE MOD_Mesh_Vars           ,ONLY:nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: VarNames(2*nVar)
CHARACTER(LEN=255)            :: FileName= 'SOLUTION_State.h5'
VarNames( 1)='Density'
VarNames( 2)='MomentumX'
VarNames( 3)='MomentumY'
VarNames( 4)='MomentumZ'
VarNames( 5)='EnergyStagnationDensity'
VarNames( 6)='VelocityX'
VarNames( 7)='Velocityy'
VarNames( 8)='VelocityZ'
VarNames( 9)='Pressure'
VarNames(10)='Temperature'
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND STDDEV TO HDF5 FILE "',TRIM(FileName),'" ... \n'
nGlobalElems=nElemsNew
CALL GenerateFileSkeleton(TRIM(FileName),'State',1,NNew,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                          'mesh_1.h5',Time_State,Time_State,withUserblock=.FALSE.,batchMode=.FALSE.,create=.TRUE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar+1,NNew,VarNames,&
                          'mesh_1.h5',Time_State,Time_State,create=.FALSE.,Dataset='Mean',batchMode=.FALSE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar+1,NNew,VarNames,&
                          'mesh_1.h5',Time_State,Time_State,create=.FALSE.,Dataset='StandardDeviation',batchMode=.FALSE.)

CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)

CALL WriteArray('Mean',5,(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/0,0,0,0,0/),.FALSE.,RealArray=UMean)
CALL WriteArray('StandardDeviation',5,(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),(/0,0,0,0,0/),.FALSE.,RealArray=SQRT(UVar))
CALL WriteArray(DataSetName='ElemData',&
					 rank=2,&
					 nValGlobal=(/nVal_ElemData(1),nVal_ElemData(2)/),&
					 nVal=(/nVal_ElemData(1),nVal_ElemData(2)/),&
					 offset=(/0,0/),&
					 collective=.FALSE.,&
					 RealArray=ElemData)
CALL WriteAttribute(File_ID,'VarNamesAdd',nVal_ElemData(1),StrArray=VarNamesAdditional)
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5


END MODULE MOD_Nisp_Output
