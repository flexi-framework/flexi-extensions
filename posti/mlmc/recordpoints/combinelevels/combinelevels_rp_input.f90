#include "flexi.h"

!===================================================================================================================================
!>
!===================================================================================================================================
MODULE MOD_CombineLevels_RP_Input
! MODULES
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitReadSumsFromHDF5
  MODULE PROCEDURE InitReadSumsFromHDF5
END INTERFACE

INTERFACE ReadSumsFromHDF5
  MODULE PROCEDURE ReadSumsFromHDF5
END INTERFACE

INTERFACE FinalizeReadSumsFromHDF5
  MODULE PROCEDURE FinalizeReadSumsFromHDF5
END INTERFACE

PUBLIC:: ReadSumsFromHDF5,InitReadSumsFromHDF5,FinalizeReadSumsFromHDF5
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Subroutine to write the sum of all deltaU and deltaU^2 RP data to HDF5 format for a specififc level 
!===================================================================================================================================
SUBROUTINE InitReadSumsFromHDF5()
! MODULES
USE MOD_Globals
USE MOD_CombineLevels_RP_Vars 
USE MOD_OutputRPVisu_HDF5
USE MOD_IO_HDF5,                 ONLY: File_ID,nDims,HSize
USE MOD_RPSetVisuVisu_Vars      ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist,OutputGroup,GroupNames
USE MOD_HDF5_Input
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,iPoint,GroupID
CHARACTER(LEN=255)               :: FileType_tmp
CHARACTER(LEN=255)               :: DataSetName
CHARACTER(LEN=255)               :: ZoneTitle
CHARACTER(LEN=255)               :: GroupName
LOGICAL                          :: DataSetFound=.FALSE.
!===================================================================================================================================
!WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" READ SUMS OF MLMC RP DATA FROM HDF5 FILE '",TRIM(FileNameSumsIn),"'..."
!WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" READ SUMS OF MLMC RP DATA FROM HDF5 FILE '",TRIM(FileNameMean),"'..."
CALL OpenDataFile(TRIM(FileNameMean),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)

!! Write dataset attributes
CALL ReadAttribute(File_ID,'File_Type'   ,1,StrScalar=FileType_tmp)
IF(TRIM(FileType_tmp) .NE. 'RP_Output') THEN
  WRITE(UNIT_stdOut,'(A)')' The Input RP File is not valid, skipping!'
  RETURN
END IF
CALL GetAttributeSize(File_ID,'VarNames',nDims,HSize)

nVal=INT(HSize(1),4)
ALLOCATE(VarNamesRP(nVal))
CALL ReadAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNamesRP)


WRITE(ZoneTitle,'(A)')'Time'
CALL DatasetExists(File_ID, ZoneTitle, DataSetFound)
CALL GetDataSize(File_ID,ZoneTitle,nDims,HSize)
nF=INT(HSize(1),4)
ALLOCATE(RP_Freq(nF))
CALL ReadArray(ZoneTitle,1,(/nF/),0,1,RealArray=RP_Freq)


! Points 
iPoint=1
GroupID=Points_GroupIDlist(iPoint)
GroupName=GroupNames(GroupID)

WRITE(DataSetName,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint
CALL DatasetExists(File_ID, DataSetName, DataSetFound)
CALL GetDataSize(File_ID,DataSetName,nDims,HSize)
nVal=INT(HSize(1),4)
nF=INT(HSize(2),4)
ALLOCATE(Utmp(nVal,nPoints,nF))
ALLOCATE(CoordinatesRP(3,nPoints))

! Close the file.
CALL CloseDataFile()

!ALLOCATE(USum(nVal/2,nPoints,nF))
!ALLOCATE(USqSum(nVal/2,nPoints,nF))
ALLOCATE(UMean(nVal,nPoints,nF))
ALLOCATE(UVariance(nVal,nPoints,nF))
UMean=0.
UVariance=0.


END SUBROUTINE InitReadSumsFromHDF5



!===================================================================================================================================
!> Subroutine to write the sum of all deltaU and deltaU^2 RP data to HDF5 format for a specififc level 
!===================================================================================================================================
SUBROUTINE ReadSumsFromHDF5()
! MODULES
USE MOD_Globals
USE MOD_CombineLevels_RP_Vars 
USE MOD_OutputRPVisu_HDF5
USE MOD_IO_HDF5,                 ONLY: File_ID,nDims,HSize
USE MOD_RPSetVisuVisu_Vars      ,ONLY:nPoints,Points_IDlist,Points_GroupIDlist,OutputGroup,GroupNames
USE MOD_HDF5_Input
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: i,iPoint,GroupID
REAL,ALLOCATABLE                 :: tmp(:,:,:)
CHARACTER(LEN=255)               :: FileType_tmp
CHARACTER(LEN=255)               :: DataSetName
CHARACTER(LEN=255)               :: ZoneTitle
CHARACTER(LEN=255)               :: GroupName
LOGICAL                          :: DataSetFound=.FALSE.
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" READ MEAN OF MLMC RP SPEC DATA FROM HDF5 FILE '",TRIM(FileNameMean),"'...\n"
CALL OpenDataFile(TRIM(FileNameMean),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)


! Points 
DO iPoint=1,nPoints
  GroupID=Points_GroupIDlist(iPoint)
  GroupName=GroupNames(GroupID)

  !values
  WRITE(DataSetName,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint
  CALL ReadArray(DataSetName,2,(/nVal,nF/),0,2,RealArray=Utmp(:,iPoint,:))
  
  !coordinates
  WRITE(DataSetName,'(A,A,I0.4,A)')TRIM(GroupName),'_Point',iPoint,'_X'
  CALL ReadArray(DataSetName,1,(/3,1/),0,1,RealArray=CoordinatesRP(:,iPoint))

END DO ! iPoint
! Close the file.
CALL CloseDataFile()

UMean(:,:,:)=UMean(:,:,:)+Utmp(:,:,:) 

WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" READ VARIANCE OF MLMC RP SPEC DATA FROM HDF5 FILE '",TRIM(FileNameVariance),"'...\n"
CALL OpenDataFile(TRIM(FileNameVariance),create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)


! Points 
DO iPoint=1,nPoints
  GroupID=Points_GroupIDlist(iPoint)
  GroupName=GroupNames(GroupID)

  !values
  WRITE(DataSetName,'(A,A,I0.4)')TRIM(GroupName),'_Point',iPoint
  CALL ReadArray(DataSetName,2,(/nVal,nF/),0,2,RealArray=Utmp(:,iPoint,:))

END DO ! iPoint
! Close the file.
CALL CloseDataFile()

UVariance(:,:,:)=UVariance(:,:,:)+Utmp(:,:,:) 

END SUBROUTINE ReadSumsFromHDF5


!===================================================================================================================================
!> Subroutine to write the sum of all deltaU and deltaU^2 RP data to HDF5 format for a specififc level 
!===================================================================================================================================
SUBROUTINE FinalizeReadSumsFromHDF5()
! MODULES
USE MOD_Globals
USE MOD_CombineLevels_RP_Vars 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

SDEALLOCATE(VarNamesRP)
SDEALLOCATE(CoordinatesRP)
SDEALLOCATE(RP_Freq)
SDEALLOCATE(Utmp)
SDEALLOCATE(UMean)
SDEALLOCATE(UVariance) 

END SUBROUTINE FinalizeReadSumsFromHDF5

END MODULE MOD_CombineLevels_RP_Input



