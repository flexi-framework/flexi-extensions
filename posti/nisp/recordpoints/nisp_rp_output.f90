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

PUBLIC:: WriteMeanAndVarianceToHDF5
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
USE MOD_Nisp_RP_Vars        ,ONLY: time, FileNameTS,FileNameFFT, UMeanTimeseries,UMeanFFT, UVarTimeseries, UVarFFT
USE MOD_Nisp_RP_Vars        ,ONLY: UTimeseries,nStochSamples,doSPL,UFFT
USE MOD_spec_Vars           ,ONLY: nSamples_spec,RPData_freq
USE MOD_ParametersVisu      ,ONLY: nVarVisu,VarNameVisu,ProjectName
USE MOD_OutputRPVisu_HDF5   ,ONLY: WriteDataToHDF5
USE MOD_HDF5_Output         ,ONLY: WriteAttribute,WriteArray,GenerateFileSkeleton
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                          :: nVarVisu_loc
CHARACTER(LEN=255),ALLOCATABLE   :: VarNames_tmp(:)
REAL,ALLOCATABLE                 :: MeanAndVarianceTS(:,:,:,:)
REAL,ALLOCATABLE                 :: MeanAndVarianceSpec(:,:,:,:)
!===================================================================================================================================
! Output spectra
!CoordNames(1)='Frequency'
IF(doSPL) THEN
  nVarVisu_loc=nVarVisu+1
ELSE
  nVarVisu_loc = nVarVisu
END IF

ALLOCATE (VarNames_tmp(nVarVisu_loc))
VarNames_tmp(1:nVarVisu)=VarNameVisu
ALLOCATE (MeanAndVarianceSpec(nVarVisu_loc,nPoints,nSamples_spec,2))
ALLOCATE (MeanAndVarianceTS(nVarVisu,nPoints,nSamples_out,2))
IF(doSPL) VarNames_tmp(nVarVisu_loc) = 'SPL'

!-----------------------------------------------------------------------------------------------------------------------------------
! TimeSeries: 
!-----------------------------------------------------------------------------------------------------------------------------------
FilenameTS=TRIM(ProjectName)//'_RP_PP.h5'
MeanAndVarianceTS(:,:,:,1)=UMeanTimeseries
MeanAndVarianceTS(:,:,:,2)=UVarTimeseries
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN, VARIANCE AND SAMPLE OF TIMESERIES TO ', FileNameTS
CALL WriteRPDataToHDF5(nSamples_out,nPoints,nVarVisu,VarNames_tmp,time,MeanAndVarianceTS,nStochSamples,UTimeseries,FileNameTS)
WRITE(UNIT_StdOut,'(132("-"))')

!-----------------------------------------------------------------------------------------------------------------------------------
! Spec: 
!-----------------------------------------------------------------------------------------------------------------------------------
FilenameFFT=TRIM(ProjectName)//'_RP_PP_spec.h5'
MeanAndVarianceSpec(:,:,:,1)=UMeanFFT
MeanAndVarianceSpec(:,:,:,2)=UVarFFT
WRITE(UNIT_StdOut,'(132("-"))')
WRITE(UNIT_stdOut,'(A,A)')' WRITING MEAN, VARIANCE AND SAMPLE OF FFT TO ', FileNameFFT
CALL WriteRPDataToHDF5(nSamples_spec,nPoints,nVarVisu_loc,VarNames_tmp,RPData_freq,MeanAndVarianceSpec,nStochSamples,UFFT,FileNameFFT)
WRITE(UNIT_StdOut,'(132("-"))')

SDEALLOCATE (VarNames_tmp)
SDEALLOCATE (MeanAndVarianceTS)
SDEALLOCATE (MeanAndVarianceSpec)

END SUBROUTINE WriteMeanAndVarianceToHDF5

!===================================================================================================================================
!> Subroutine to write record point data to HDF5 file
!===================================================================================================================================
SUBROUTINE WriteRPDataToHDF5(nSamples,nRP,nVal,VarNames,Time,MeanAndVariance,nStSamples,Value,FileString)
! MODULES
USE MOD_Globals
USE HDF5
USE MOD_IO_HDF5
USE MOD_HDF5_Output
USE MOD_ParametersVisu     ,ONLY: ProjectName
USE MOD_ParametersVisu     ,ONLY: Line_LocalCoords,Plane_LocalCoords
USE MOD_ParametersVisu     ,ONLY: OutputPlanes,OutputLines,OutputPoints
USE MOD_RPSetVisuVisu_Vars ,ONLY: GroupNames
USE MOD_RPSetVisuVisu_Vars ,ONLY: OutputGroup
USE MOD_RPSetVisuVisu_Vars ,ONLY: nPoints,Points_IDlist,Points_GroupIDlist
USE MOD_RPSetVisuVisu_Vars ,ONLY: nLines,Lines,tLine
USE MOD_RPSetVisuVisu_Vars ,ONLY: nPlanes,Planes,tPlane
USE MOD_RPSetVisuVisu_Vars ,ONLY: xF_RP
USE MOD_OutputRPVisu_Vars  ,ONLY: nCoords,CoordNames
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: nSamples                      !< Number of Samples
INTEGER,INTENT(IN)            :: nRP                           !< Number of RP to be visualized
INTEGER,INTENT(IN)            :: nVal                          !< Number of nodal output variables
CHARACTER(LEN=255),INTENT(IN) :: VarNames(nVal)                !< Names of all variables that will be written out
REAL,INTENT(IN)               :: MeanAndVariance(1:nVal,nRP,nSamples,2)    !< Statevector
INTEGER,INTENT(IN)            :: nStSamples                    !< Number of nodal output variables
REAL,INTENT(IN)               :: Value(1:nVal,nRP,nSamples,nStSamples)   !< Statevector
REAL,INTENT(IN)               :: Time(nSamples)                !< Time
CHARACTER(LEN=*),INTENT(IN)   :: FileString                    !< Output file name
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER              :: iVar
INTEGER              :: nCoords_loc
INTEGER              :: iPoint,iLine,iPlane,i,j
INTEGER              :: GroupID
INTEGER              :: size_offsetdim,offsetVar
CHARACTER(LEN=255)   :: ZoneTitle
CHARACTER(LEN=255)   :: ZoneTitleMean
CHARACTER(LEN=255)   :: GroupName
CHARACTER(LEN=255)   :: CoordNames_loc(nCoords-1)
REAL                 :: PointDataMeanAndVar(1:nVal,nSamples,2)
REAL                 :: PointData(1:nVal,nSamples,nStSamples)
TYPE(tLine),POINTER  :: Line
TYPE(tPlane),POINTER :: Plane
REAL,ALLOCATABLE     :: LineDataMeanAndVar(:,:,:,:)
REAL,ALLOCATABLE     :: PlaneDataMeanAndVar(:,:,:,:,:)
REAL,ALLOCATABLE     :: LineData(:,:,:,:)
REAL,ALLOCATABLE     :: PlaneData(:,:,:,:,:)
REAL,ALLOCATABLE     :: PlaneCoord(:,:,:)
!===================================================================================================================================
WRITE(UNIT_stdOut,'(A,A,A)',ADVANCE='NO')" WRITE RP DATA TO HDF5 FILE '",TRIM(FileString),"'..."
CALL OpenDataFile(Filestring,create=.TRUE.,single=.TRUE.,readOnly=.FALSE.)

! Write dataset attributes
CALL WriteAttribute(File_ID,'File_Type'   ,1,StrScalar=(/'RP_Output'/))
CALL WriteAttribute(File_ID,'ProjectName' ,1,StrScalar=(/ProjectName/))
CALL WriteAttribute(File_ID,'Time'        ,1,RealScalar=Time(nSamples))
CALL WriteAttribute(File_ID,'VarNames'    ,nVal,StrArray=VarNames)
nCoords_loc=nCoords-1
CoordNames_loc=CoordNames(2:nCoords)
CALL WriteAttribute(File_ID,'CoordNames',nCoords_loc,StrArray=CoordNames_loc)

!time (global for all points)
IF(nSamples.GT.1)THEN
  WRITE(ZoneTitle,'(A)')'Time'
  offsetVar=0
  CALL WriteArray(TRIM(ZoneTitle),1,(/nSamples/),(/nSamples/),(/0/),.FALSE.,RealArray=Time)
END IF

! Points
IF(OutputPoints)THEN
  DO iPoint=1,nPoints
    GroupID=Points_GroupIDlist(iPoint)
    IF(.NOT.OutputGroup(GroupID)) CYCLE
    GroupName=GroupNames(GroupID)

    ! coordinates
    ZoneTitle(1:255)=' '
    WRITE(ZoneTitle,'(A,A,I0.4,A)')TRIM(GroupName),'_Point_',iPoint,'_X'
    CALL WriteArray(TRIM(ZoneTitle),1,(/3/),(/3/),(/0/),.FALSE.,RealArray=xF_RP(:,Points_IDlist(iPoint)))

    !Mean and Variance
    WRITE(ZoneTitleMean,'(A,A,I0.4,A)')TRIM(GroupName),'_Point_',iPoint,'_MeanAndVariance'
    WRITE(ZoneTitle,'(A,A,I0.4,A)')TRIM(GroupName),'_Point_',iPoint,'_Samples'
    size_offsetdim=nVal
    offsetVar=0
    PointDataMeanAndVar(:,:,:)=MeanAndVariance(:,Points_IDlist(iPoint),:,:)
    PointData(:,:,:)=Value(:,Points_IDlist(iPoint),:,:)
    DO iVar=1,nVal
      CALL WriteArray(TRIM(ZoneTitleMean),3,(/nVal,nSamples,2/),(/1,nSamples,2/),(/offsetVar,0,0/),&
                      .FALSE.,RealArray=PointDataMeanAndVar(iVar,:,:))
      CALL WriteArray(TRIM(ZoneTitle),3,(/nVal,nSamples,nStSamples/),(/1,nSamples,nStSamples/),(/offsetVar,0,0/),&
                      .FALSE.,RealArray=PointData(iVar,:,:))
      offsetVar=offsetVar+1
    END DO       ! iVar
  END DO ! iPoint
END IF!OutputPoint

! Lines
IF(OutputLines)THEN
  DO iLine=1,nLines
    Line=>Lines(iLine)
    IF(.NOT.OutputGroup(Line%GroupID)) CYCLE
    GroupName=GroupNames(Line%GroupID)

    !coordinates
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name),'_X'
    size_offsetdim=3
    OffsetVar=0
    ! local coord if required
    IF(Line_LocalCoords) THEN
      size_offsetdim=size_offsetdim+1
      CALL WriteArray(TRIM(ZoneTitle),2,(/size_offsetDim,Line%nRP/),(/1,Line%nRP/),(/OffsetVar,0/),.FALSE.,RealArray=Line%LocalCoord(:))
      OffsetVar=1
    END IF!(Line_LocalCoords)
    CALL WriteArray(  TRIM(ZoneTitle),2,(/size_offsetDim,Line%nRP/),(/3,Line%nRP/),(/OffsetVar,0/),.FALSE.,RealArray=xF_RP(:,Line%IDlist(:)))

    !values
    WRITE(ZoneTitleMean,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name),'_MeanAndVariance'
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Line%Name),'_Samples'
    size_offsetdim=nVal
    OffsetVar=0
    ALLOCATE(LineDataMeanAndVar(1:nVal,nSamples,Line%nRP,2))
    ALLOCATE(LineData(1:nVal,nSamples,Line%nRP,nStSamples))
    DO iPoint=1,Line%nRP
      LineDataMeanAndVar(:,:,iPoint,:)=MeanAndVariance(:,Line%IDlist(iPoint),:,:)
      LineData(:,:,iPoint,:)=Value(:,Line%IDlist(iPoint),:,:)
    END DO ! iPoint
    CALL WriteArray(TRIM(ZoneTitle),4,(/nVal,nSamples,Line%nRP,2/),(/nVal,nSamples,Line%nRP,2/),(/0,0,0,0/),.FALSE.,RealArray=LineDataMeanAndVar)
    CALL WriteArray(TRIM(ZoneTitle),4,(/nVal,nSamples,Line%nRP,nStSamples/),(/nVal,nSamples,Line%nRP,nStSamples/),(/0,0,0,0/),.FALSE.,RealArray=LineData)
    DEALLOCATE(LineDataMeanAndVar)
    DEALLOCATE(LineData)
  END DO ! iLine
END IF !OutputLines

! Planes
IF(OutputPlanes)THEN
  DO iPlane=1,nPlanes
    Plane=>Planes(iPlane)
    IF(.NOT.OutputGroup(Plane%GroupID)) CYCLE
    GroupName=GroupNames(Plane%GroupID)
    ZoneTitle(1:255)=' '


    !coordinates
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_X'
    size_offsetdim=3
    OffsetVar=0
    ! local coord if required
    IF(Plane_LocalCoords.AND.ALLOCATED(Plane%LocalCoord)) THEN
      size_offsetdim=size_offsetdim+2
      CALL WriteArray(TRIM(ZoneTitle),3,(/size_offsetdim,Plane%nRP(1),Plane%nRP(2)/),(/2,Plane%nRP(1),Plane%nRP(2)/),(/OffsetVar,0,0/),.FALSE.,RealArray=Plane%LocalCoord(:,:,:))
      OffsetVar=2
    END IF!(Line_LocalCoords)
    ALLOCATE(PlaneCoord(3,Plane%nRP(1),Plane%nRP(2)))
    ! global xyz coordinates
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneCoord(:,i,j)=xF_RP(:,Plane%IDlist(i,j))
      END DO ! i
    END DO ! j
    CALL WriteArray(TRIM(ZoneTitle),3,(/size_offsetdim,Plane%nRP(1),Plane%nRP(2)/),(/3,Plane%nRP(1),Plane%nRP(2)/),&
                                 (/OffsetVar,0,0/),.FALSE.,RealArray=PlaneCoord)
    DEALLOCATE(PlaneCoord)

    !values
    WRITE(ZoneTitleMean,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_MeanAndVariance'
    WRITE(ZoneTitle,'(A,A,A,A)')TRIM(GroupName),'_',TRIM(Plane%Name),'_Samples'
    size_offsetdim=nVal
    OffsetVar=0
    ALLOCATE(PlaneDataMeanAndVar(1:nVal,nSamples,Plane%nRP(1),Plane%nRP(2),2))
    ALLOCATE(PlaneData(1:nVal,nSamples,Plane%nRP(1),Plane%nRP(2),nStSamples))
    DO j=1,Plane%nRP(2)
      DO i=1,Plane%nRP(1)
        PlaneDataMeanAndVar(:,:,i,j,:)=MeanAndVariance(:,Plane%IDlist(i,j),:,:)
        PlaneData(:,:,i,j,:)=Value(:,Plane%IDlist(i,j),:,:)
      END DO ! i
    END DO ! j
    CALL WriteArray(TRIM(ZoneTitleMean),5,(/nVal,nSamples,Plane%nRP(1),Plane%nRP(2),2/),(/nVal,nSamples,Plane%nRP(1),Plane%nRP(2),2/),&
                    (/0,0,0,0,0/),.FALSE.,RealArray=PlaneDataMeanAndVar)
    CALL WriteArray(TRIM(ZoneTitle),5,(/nVal,nSamples,Plane%nRP(1),Plane%nRP(2),nStSamples/),(/nVal,nSamples,Plane%nRP(1),Plane%nRP(2),nStSamples/),&
                    (/0,0,0,0,0/),.FALSE.,RealArray=PlaneDataMeanAndVar)
    DEALLOCATE(PlaneDataMeanAndVar)
    DEALLOCATE(PlaneData)
  END DO ! iPlane
END IF!OutputPlanes

! Close the file.
CALL CloseDataFile()

WRITE(UNIT_stdOut,'(A)',ADVANCE='YES')"DONE"
END SUBROUTINE WriteRPDataToHDF5

END MODULE MOD_Nisp_RP_Output
