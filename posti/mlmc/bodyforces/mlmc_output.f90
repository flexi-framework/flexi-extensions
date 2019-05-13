#include "flexi.h"

MODULE MOD_MLMC_Output
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

INTERFACE WriteBodyForcesSumsToHDF5
  MODULE PROCEDURE WriteBodyForcesSumsToHDF5
END INTERFACE

INTERFACE WriteMeanAndVarianceToHDF5
  MODULE PROCEDURE WriteMeanAndVarianceToHDF5
END INTERFACE

PUBLIC::WriteSumsToHDF5
PUBLIC::WriteBodyForcesSumsToHDF5
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
USE MOD_MLMC_Vars
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
CHARACTER(LEN=255)   :: DataSetNames(5)
INTEGER              :: i
TYPE tArrayPtr
  REAL,POINTER       :: ptr(:,:,:,:,:)
END TYPE
TYPE(tArrayPtr)                       :: DataSetPtrs(5)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE Sums TO HDF5 FILE "',TRIM(FileNameSums),'" ...'

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

CALL GenerateFileSkeleton(TRIM(FileNameSums),'EstSigSum',nVarTotal,NNew,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                          MeshFileNew,Time_State,Time_State,withUserblock=.FALSE.,batchMode=.FALSE.,create=.TRUE.)

DataSetNames(1)='UFineSum'
DataSetNames(2)='UCoarseSum'
DataSetNames(3)='UFineSqSum'
DataSetNames(4)='UCoarseSqSum'
DataSetNames(5)='DUSqSum'

DO i=1,5
  CALL GenerateFileSkeleton(TRIM(FileNameSums),'EstSigSum',nVarTotal,NNew,VarNames,&
                            MeshFileNew,Time_State,Time_State,create=.FALSE.,Dataset=DataSetNames(i),batchMode=.FALSE.)
END DO

DataSetPtrs(1)%ptr=>UFineSum
DataSetPtrs(2)%ptr=>UCoarseSum
DataSetPtrs(3)%ptr=>UFineSqSum
DataSetPtrs(4)%ptr=>UCoarseSqSum
DataSetPtrs(5)%ptr=>DUSqSum

CALL OpenDataFile(TRIM(FileNameSums),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
DO i=1,5
  CALL WriteArray(DataSetName=DataSetNames(i),&
                  rank=5,&
                  nValGlobal=(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
                  nVal=(/nVarTotal,NNew+1,NNew+1,NNew+1,nElemsNew/),&
                  offset=(/0,0,0,0,0/),&
                  collective=.FALSE.,&
                  RealArray=DataSetPtrs(i)%ptr)
END DO

CALL WriteAttribute(File_ID,'nSamples',1,IntScalar=NEnd)

CALL WriteAttribute(File_ID,'SigmaSq',1,RealScalar=SigmaSq)
CALL WriteAttribute(File_ID,'SigmaSqFine',1,RealScalar=SigmaSqFine)
CALL WriteAttribute(File_ID,'Bias',1,RealScalar=Bias)

CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteSumsToHDF5


SUBROUTINE WriteBodyForcesSumsToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_MLMC_Vars
USE MOD_HDF5_Output,                ONLY:WriteAttribute, WriteArray
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)   :: DataSetNames(5)
INTEGER              :: i
TYPE tArrayPtr
  REAL,POINTER       :: ptr(:)
END TYPE
TYPE(tArrayPtr)      :: DataSetPtrs(5)
!===================================================================================================================================
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE BODYFORCES Sums TO HDF5 FILE "',TRIM(FileNameSumsBF),'" ...'
DataSetNames(1)='BodyForcesFineSum'
DataSetNames(2)='BodyForcesCoarseSum'
DataSetNames(3)='BodyForcesFineSqSum'
DataSetNames(4)='BodyForcesCoarseSqSum'
DataSetNames(5)='DBodyForcesSqSum'
DataSetPtrs(1)%ptr=>BodyForcesFineSum
DataSetPtrs(2)%ptr=>BodyForcesCoarseSum
DataSetPtrs(3)%ptr=>BodyForcesFineSqSum
DataSetPtrs(4)%ptr=>BodyForcesCoarseSqSum
DataSetPtrs(5)%ptr=>DBodyForcesSqSum

CALL OpenDataFile(FileNameSumsBF,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

DO i=1,5
  CALL WriteArray(DataSetName=DataSetNames(i),&
                  rank=1,&
                  nValGlobal=(/1/),&
                  nVal=(/1/),&
                  offset=(/0/),&
                  collective=.FALSE.,&
                  RealArray=DataSetPtrs(i)%ptr)
END DO

CALL WriteAttribute(File_ID,'nSamples',1,IntScalar=NEnd)

CALL WriteAttribute(File_ID,'SigmaSq',1,RealScalar=SigmaSq)
CALL WriteAttribute(File_ID,'SigmaSqFine',1,RealScalar=SigmaSqFine)
CALL WriteAttribute(File_ID,'Bias',1,RealScalar=Bias)

CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteBodyForcesSumsToHDF5

SUBROUTINE WriteMeanAndVarianceToHDF5()
!===================================================================================================================================
! Subroutine to write the time averaged solution U to HDF5 format
!===================================================================================================================================
! MODULES
USE MOD_Globals
USE MOD_IO_HDF5
USE MOD_HDF5_Output,ONLY:WriteAttribute,WriteArray,GenerateFileSkeleton
USE MOD_MLMC_Vars
USE MOD_Mesh_Vars,ONLY:nGlobalElems
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)  :: FileName
INTEGER             :: NLoc
!===================================================================================================================================
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
!-----------------------------------------------------------------------------------------------------------------------------------
FileName = 'SOLUTION_State.h5'
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' WRITE MEAN AND STDDEV TO HDF5 FILE "',TRIM(FileName),'" ... \n'

NLoc=nVal(2)-1
nGlobalElems=nVal(5)
CALL GenerateFileSkeleton(TRIM(FileName),'EstSigSum',1,NLoc,(/'DUMMY_DO_NOT_VISUALIZE'/),&
                          MeshFile_Sums,Time_Sums,Time_Sums,withUserblock=.FALSE.,batchMode=.FALSE.,create=.TRUE.)
CALL GenerateFileSkeleton(TRIM(FileName),'EstSigSum',nVarTotal,NLoc,VarNames,&
                          MeshFile_Sums,Time_Sums,Time_Sums,create=.FALSE.,Dataset='Mean',batchMode=.FALSE.)
CALL GenerateFileSkeleton(TRIM(FileName),'EstSigSum',nVarTotal,NLoc,VarNames,&
                          MeshFile_Sums,Time_Sums,Time_Sums,create=.FALSE.,Dataset='StandardDeviation',batchMode=.FALSE.)
CALL OpenDataFile(TRIM(FileName),create=.FALSE.,single=.TRUE.,readOnly=.FALSE.)
CALL WriteArray('Mean',5,nVal,nVal,(/0,0,0,0,0/),.FALSE.,RealArray=Mean)
CALL WriteArray('StandardDeviation',5,nVal,nVal,(/0,0,0,0,0/),.FALSE.,RealArray=StdDev)
CALL CloseDataFile()

SWRITE(UNIT_stdOut,'(a)',ADVANCE='YES')'DONE'
END SUBROUTINE WriteMeanAndVarianceToHDF5

END MODULE MOD_MLMC_Output
