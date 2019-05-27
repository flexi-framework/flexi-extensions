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
USE MOD_Mesh_Vars           ,ONLY: nGlobalElems
USE MOD_Nisp_Vars
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
CHARACTER(LEN=255)            :: FileName
Vars( 1)='Density'
Vars( 2)='MomentumX'
Vars( 3)='MomentumY'
Vars( 4)='MomentumZ'
Vars( 5)='EnergyStagnationDensity'
Vars( 6)='VelocityX'
Vars( 7)='Velocityy'
Vars( 8)='VelocityZ'
Vars( 9)='Pressure'
Vars(10)='Temperature'
!===================================================================================================================================
!-----------------------------------------------------------------------------------------------------------------------------------
!hack: GenerateFileSkeleton needs PP_N and nGlobalElems
PP_N = NNew
nGlobalElems = nElemsNew
! FileNameMean=TRIM(TIMESTAMP('Mean',Time_State))//'.h5'
FileName = 'SOLUTION_State.h5'
!-----------------------------------------------------------------------------------------------------------------------------------
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND Variace TO HDF5 FILE "',TRIM(FileName),'" ... \n'

CALL GenerateFileSkeleton(TRIM(FileName),'State',1,NNew,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                          'mesh.h5',Time_State,Time_State,withUserblock=.FALSE.,batchMode=.FALSE.,create=.TRUE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar,NNew,Vars,&
                          'mesh.h5',Time_State,Time_State,create=.FALSE.,Dataset='Mean',batchMode=.FALSE.)
CALL GenerateFileSkeleton(TRIM(FileName),'State',2*nVar,NNew,Vars,&
                          'mesh.h5',Time_State,Time_State,create=.FALSE.,Dataset='Variance',batchMode=.FALSE.)
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
CALL WriteArray('Mean',5,nVal=(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),nValGlobal=(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),offset=(/0,0,0,0,0/),collective=.FALSE.,RealArray=UMean)
CALL WriteArray('Variance',5,nVal=(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),nValGlobal=(/2*nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),offset=(/0,0,0,0,0/),collective=.FALSE.,RealArray=UVar)
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5


END MODULE MOD_Nisp_Output
