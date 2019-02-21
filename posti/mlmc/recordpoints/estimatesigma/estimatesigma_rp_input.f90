#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_EstimateSigma_RP_Input
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE ReadSumsFromHDF5
  MODULE PROCEDURE ReadSumsFromHDF5
END INTERFACE

PUBLIC:: ReadSumsFromHDF5
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Subroutine to write the sum of all deltaU and deltaU^2 RP data to HDF5 format for a specififc level 
!===================================================================================================================================
SUBROUTINE ReadSumsFromHDF5()
! MODULES
USE MOD_Globals
!USE MOD_RPData_Vars         ,ONLY: RPTime
!USE MOD_ParametersVisu      ,ONLY: OutputFormat,thirdOct,ProjectName
!USE MOD_OutputRPVisu_Vars   ,ONLY: nSamples_out,RPData_out,RPDataTimeAvg_out,CoordNames 
!USE MOD_OutputRPVisu_HDF5
!USE MOD_spec_Vars           ,ONLY: nSamples_spec,RPData_freq,RPData_spec
!USE MOD_spec_Vars           ,ONLY: nSamples_Oct,RPData_freqOct,RPData_Oct
!USE MOD_ParametersVisu      ,ONLY: nVarVisu,VarNameVisu
!USE MOD_ParametersVisu      ,ONLY: OutputTimeAverage,OutputTimeData,doSpec,doFluctuations
!USE MOD_ParametersVisu      ,ONLY: Plane_doBLProps
!USE MOD_RPSetVisuVisu_Vars  ,ONLY: nRP_global
USE MOD_EstimateSigma_RP_Vars 
USE MOD_spec_Vars             ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu        ,ONLY: nVarVisu,VarNameVisu
USE MOD_OutputRPVisu_HDF5
USE MOD_IO_HDF5,                 ONLY: File_ID,nDims,HSize
USE MOD_RPSetVisuVisu_Vars      ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist,OutputGroup,GroupNames
USE MOD_HDF5_Input
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,nVal,iPoint,GroupID,nF
REAL,ALLOCATABLE                 :: tmp(:,:,:)
REAL,ALLOCATABLE                 :: SumsOld(:,:)
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:)
CHARACTER(LEN=255)               :: FileType_tmp
CHARACTER(LEN=255)               :: DataSetName
CHARACTER(LEN=255)               :: GroupName
LOGICAL                          :: DataSetFound=.FALSE.
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" READ SUMS OF MLMC RP DATA FROM HDF5 FILE '",TRIM(FileNameSumsIn),"'..."
CALL OpenDataFile(FileNameSumsIn,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

!! Write dataset attributes
CALL ReadAttribute(File_ID,'File_Type'   ,1,StrScalar=FileType_tmp)
IF(TRIM(FileType_tmp) .NE. 'RP_Output') THEN
  WRITE(UNIT_stdOut,'(A)')' The Input RP File is not valid, skipping!'
  RETURN
END IF
CALL GetAttributeSize(File_ID,'VarNames',nDims,HSize)

nVal=INT(HSize(1),4)
ALLOCATE(VarNames_tmp(nVal))
CALL ReadAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNames_tmp)


! Points 
DO iPoint=1,nPoints
  GroupID=Points_GroupIDlist(iPoint)
  IF(.NOT.OutputGroup(GroupID)) CYCLE
  GroupName=GroupNames(GroupID)


  !values
  WRITE(DataSetName,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint
  CALL DatasetExists(File_ID, DataSetName, DataSetFound)
  CALL GetDataSize(File_ID,DataSetName,nDims,HSize)
  
  nVal=INT(HSize(1),4)
  nF=INT(HSize(2),4)
  ALLOCATE(SumsOld(nVal,nF))
  CALL ReadArray(DataSetName,2,(/nVal,nF/),0,2,RealArray=SumsOld)
  
  !USum(:,iPoint,:)=SumsOld(1:(nVal/2),:) 
  !USqSum(:,iPoint,:)=SumsOld((nVal/2)+1:nVal,:) 
  nVal=nVal/5
  UFineSum(:,iPoint,:)     = SumsOld(1       :  nVal,:)
  UCoarseSum(:,iPoint,:)   = SumsOld(  nVal+1:2*nVal,:)
  UFineSqSum(:,iPoint,:)   = SumsOld(2*nVal+1:3*nVal,:)
  UCoarseSqSum(:,iPoint,:) = SumsOld(3*nVal+1:4*nVal,:)
  DUSqSum(:,iPoint,:)      = SumsOld(4*nVal+1:      ,:)

  SDEALLOCATE(SumsOld)
END DO ! iPoint

! Close the file.
CALL CloseDataFile()

END SUBROUTINE ReadSumsFromHDF5


END MODULE MOD_EstimateSigma_RP_Input



