#include "flexi.h"

MODULE MOD_EstimateSigma_Output
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

INTERFACE WriteSumsToHDF5
  MODULE PROCEDURE WriteSumsToHDF5
END INTERFACE

INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

PUBLIC::WriteSumsToHDF5
PUBLIC::WriteMeanAndVarianceToHDF5
!===================================================================================================================================
CONTAINS



SUBROUTINE WriteSumsToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_EstimateSigma_Vars
USE MOD_HDF5_Output,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_SwapMesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,ALLOCATABLE     :: tmp(:,:,:,:,:)
CHARACTER(LEN=255)   :: VarNames(5*nVarTotal)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE Sums TO HDF5 FILE "',TRIM(FileNameSums),'" ...'

VarNames( 1)=                'Density_F'
VarNames( 2)=              'MomentumX_F'
VarNames( 3)=              'MomentumY_F'
VarNames( 4)=              'MomentumZ_F'
VarNames( 5)='EnergyStagnationDensity_F'
VarNames( 6)=              'VelocityX_F'
VarNames( 7)=              'Velocityy_F'
VarNames( 8)=              'VelocityZ_F'
VarNames( 9)=               'Pressure_F'
VarNames(10)=            'Temperature_F'

VarNames(11)=                'Density_C'
VarNames(12)=              'MomentumX_C'
VarNames(13)=              'MomentumY_C'
VarNames(14)=              'MomentumZ_C'
VarNames(15)='EnergyStagnationDensity_C'
VarNames(16)=              'VelocityX_C'
VarNames(17)=              'Velocityy_C'
VarNames(18)=              'VelocityZ_C'
VarNames(19)=               'Pressure_C'
VarNames(20)=            'Temperature_C'

VarNames(21)=                'Density_FSq'
VarNames(22)=              'MomentumX_FSq'
VarNames(23)=              'MomentumY_FSq'
VarNames(24)=              'MomentumZ_FSq'
VarNames(25)='EnergyStagnationDensity_FSq'
VarNames(26)=              'VelocityX_FSq'
VarNames(27)=              'Velocityy_FSq'
VarNames(28)=              'VelocityZ_FSq'
VarNames(29)=               'Pressure_FSq'
VarNames(30)=            'Temperature_FSq'

VarNames(31)=                'Density_CSq'
VarNames(32)=              'MomentumX_CSq'
VarNames(33)=              'MomentumY_CSq'
VarNames(34)=              'MomentumZ_CSq'
VarNames(35)='EnergyStagnationDensity_CSq'
VarNames(36)=              'VelocityX_CSq'
VarNames(37)=              'Velocityy_CSq'
VarNames(38)=              'VelocityZ_CSq'
VarNames(39)=               'Pressure_CSq'
VarNames(40)=            'Temperature_CSq'

VarNames(41)=                'Density_DSq'
VarNames(42)=              'MomentumX_DSq'
VarNames(43)=              'MomentumY_DSq'
VarNames(44)=              'MomentumZ_DSq'
VarNames(45)='EnergyStagnationDensity_DSq'
VarNames(46)=              'VelocityX_DSq'
VarNames(47)=              'Velocityy_DSq'
VarNames(48)=              'VelocityZ_DSq'
VarNames(49)=               'Pressure_DSq'
VarNames(50)=            'Temperature_DSq'


ALLOCATE(tmp(5*nVarTotal,0:NNew,0:NNew,0:NNew,nElemsNew))
!-----------------------------------------------------------------------------------------------------------------------------------
tmp(1            :  nVarTotal,:,:,:,:) = UFineSum
tmp(  nVarTotal+1:2*nVarTotal,:,:,:,:) = UCoarseSum
tmp(2*nVarTotal+1:3*nVarTotal,:,:,:,:) = UFineSqSum
tmp(3*nVarTotal+1:4*nVarTotal,:,:,:,:) = UCoarseSqSum
tmp(4*nVarTotal+1:           ,:,:,:,:) = DUSqSum

CALL GenerateFileSkeleton(TRIM(FileNameSums),&
                          'EstSigSum',&
                          5*nVarTotal,&
                          NNew,&
                          VarNames,&
                          MeshFileNew,&
                          Time_State,&
                          Time_State,&
                          withUserblock=.FALSE.,&
                          batchMode=.FALSE.)
CALL OpenDataFile(TRIM(FileNameSums),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
!SUBROUTINE OpenDataFile(FileString,create,single,readOnly,communicatorOpt,userblockSize)
CALL WriteArray(DataSetName='DG_Solution',&
                rank=5,&
                nValGlobal=(/5*nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
                nVal=(/5*nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
                offset=(/0,0,0,0,0/),&
                collective=.FALSE.,&
                RealArray=tmp)
!SUBROUTINE WriteArray(DataSetName,rank,nValGlobal,nVal,offset,&
                            !collective,resizeDim,chunkSize,&
                            !RealArray,IntArray,StrArray)
CALL WriteAttribute(File_ID,'NSamples',1,IntScalar=NEnd)

CALL WriteAttribute(File_ID,'SigmaSq',1,RealScalar=SigmaSq)
CALL WriteAttribute(File_ID,'SigmaSqFine',1,RealScalar=SigmaSqFine)
CALL WriteAttribute(File_ID,'Bias',1,RealScalar=Bias)

CALL CloseDataFile()

SDEALLOCATE(tmp)

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteSumsToHDF5



SUBROUTINE WriteMeanAndVarianceToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Output,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_EstimateSigma_Vars
USE MOD_SwapMesh_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)            :: VarNames(nVarTotal)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND VARIANCE TO HDF5 FILE "',TRIM(FileNameMean),'" ...'

!-----------------------------------------------------------------------------------------------------------------------------------
! MEAN:

VarNames(1)='Density'
VarNames(2)='MomentumX'
VarNames(3)='MomentumY'
VarNames(4)='MomentumZ'
VarNames(5)='EnergyStagnationDensity'
VarNames(6)='VelocityX'
VarNames(7)='VelocityY'
VarNames(8)='VelocityZ'
VarNames(9)='Pressure'
VarNames(10)='Temperature'

CALL GenerateFileSkeleton(TRIM(FileNameMean),&
                          'State',&
                          nVarTotal,&
                          NNew,&
                          VarNames,&
                          MeshFileNew,&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)

! Open HDF5 file for write
CALL OpenDataFile(TRIM(FileNameMean),.FALSE.,.TRUE.,readOnly=.FALSE.)
! Write DG solution ----------------------------------------------------------------------------------------------------------------
CALL WriteArray('DG_Solution',5,(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    (/0,0,0,0,0/),.FALSE.,RealArray=snSamples*(UFineSum-UCoarseSum))

! Write properties -----------------------------------------------------------------------------------------------------------------

CALL CloseDataFile()



!-----------------------------------------------------------------------------------------------------------------------------------
! VARIANCE:

CALL GenerateFileSkeleton(TRIM(FileNameVariance),&
                          'State',&
                          nVarTotal,&
                          NNew,&
                          VarNames,&
                          MeshFileNew,&
                          Time_State,&
                          Time_State,&
                          withUserblock=.TRUE.)

! Open HDF5 file for write
CALL OpenDataFile(TRIM(FileNameVariance),.FALSE.,.TRUE.,readOnly=.FALSE.)
! Write DG solution ----------------------------------------------------------------------------------------------------------------
CALL WriteArray('DG_Solution',5,(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
    (/0,0,0,0,0/),.FALSE.,RealArray=SigmaSqField)

! Write properties -----------------------------------------------------------------------------------------------------------------

CALL CloseDataFile()


SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5



END MODULE MOD_EstimateSigma_Output
