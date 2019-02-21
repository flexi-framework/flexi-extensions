#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files 
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM TimeAvg
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_StringTools,ONLY:STRICMP
USE MOD_IO_HDF5
USE MOD_HDF5_Input, ONLY:ReadAttribute,ReadArray,GetDataSize
USE MOD_MPI,        ONLY:InitMPI
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
INTEGER,PARAMETER                    :: maxDim=16  !< maximum number of permitted array dimension 
TYPE tFileSet
  INTEGER                            :: nDataSets  !< number of datasets
  INTEGER(KIND=8)                    :: totalsize  !< 1D size of all arrays
  INTEGER,ALLOCATABLE                :: nDims(:)   !< number of dimensions per dataset
  INTEGER,ALLOCATABLE                :: nVal(:,:)  !< number of entries per dataset
  CHARACTER(LEN=255),ALLOCATABLE     :: DatasetNames(:) !< names of the datasets
  CHARACTER(LEN=255)                 :: FileType,MeshFile,NodeType,ProjectName
  REAL                               :: time       !< time
END TYPE
TYPE(tFileSet)                       :: ref

! LOCAL VARIABLES
CHARACTER(LEN=255)                   :: InputFile,dir
CHARACTER(LEN=255)                   :: level_string
REAL,ALLOCATABLE                     :: USum(:),ULoc(:)
INTEGER                              :: iLevel
INTEGER                              :: locsize(16),i,n
INTEGER                              :: offset,startInd,endInd
LOGICAL                              :: ex
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL DefineParametersIO_HDF5()
CALL InitMPI()
CALL InitIOHDF5()
CALL ParseCommandlineArguments()
IF ((doPrintHelp.GT.0) .OR. (nArgs.LT.1)) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: ./mlmc_combinelevels FILENAME (e.g. mean_state.h5)')
END IF

CALL GetParams('level_1/'//TRIM(Args(1)),ref)
ALLOCATE(USum(ref%totalsize))
ALLOCATE(ULoc(ref%totalsize))
USum=0.

! Start the averaging
WRITE(UNIT_stdOut,'(132("="))')

DO iLevel=1,100
  WRITE(level_string,"(I0)") iLevel
  InputFile='level_'//TRIM(level_string)//'/'//TRIM(Args(1))
  INQUIRE(FILE=TRIM(InputFile), EXIST=ex )
  IF (.NOT.ex) EXIT
  SWRITE(UNIT_stdOut,'(A,A,A)') "Processing level ",level_string,"..."

  CALL OpenDataFile(InputFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
  offset=0
  DO i=1,ref%nDataSets
    n=ref%nDims(i)
    locsize(1:n)=ref%nVal(1:n,i)
    startind=offset+1
    endind  =offset+PRODUCT(locsize(1:n))
    CALL ReadArray(TRIM(ref%DatasetNames(i)),n,locsize(1:n),0,n,&
                   RealArray=Uloc(startInd:endInd))
    offset=endInd
  END DO
  CALL CloseDataFile()

  USum = USum + Uloc

END DO
InputFile='level_1/'//TRIM(Args(1))
CALL WriteTimeAverageByCopy(InputFile,Args(1),ref,USum)
 
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') "Merging DONE!"
SWRITE(UNIT_stdOut,'(132("="))')

SDEALLOCATE(USum)
SDEALLOCATE(ULoc)
CONTAINS



!===================================================================================================================================
!> Retrieves relevant header and dateset parameters from Flexi files and stores them in a type
!===================================================================================================================================
SUBROUTINE GetParams(filename,f)
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename !< input filename
TYPE(tFileSet),INTENT(OUT)  :: f        !< type with infos to be filled
!===================================================================================================================================
CALL OpenDataFile(filename,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDatasetNamesInGroup("/",f%DatasetNames)
f%nDataSets=SIZE(f%DatasetNames)
ALLOCATE(f%nDims(f%nDataSets))
ALLOCATE(f%nVal(maxDim,f%nDataSets))
f%nVal=0
f%totalsize=0
DO i=1,f%nDataSets
  CALL GetDataSize(File_ID,TRIM(f%DatasetNames(i)),f%nDims(i),HSize)
  CHECKSAFEINT(MAXVAL(HSize),4)
  CHECKSAFEINT(MINVAL(HSize),4)
  f%nVal(1:f%nDims(i),i)=INT(HSize)
  DEALLOCATE(HSize)
  f%totalsize=f%totalsize+PRODUCT(f%nVal(1:f%nDims(i),i))
END DO
! Get default parameters
CALL ReadAttribute(File_ID,'File_Type',1,StrScalar =f%FileType)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar =f%MeshFile)
CALL ReadAttribute(File_ID,'NodeType',1,StrScalar =f%NodeType)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar =f%ProjectName)
CALL ReadAttribute(File_ID,'Time'    ,1,RealScalar=f%Time)
CALL CloseDataFile()
END SUBROUTINE


!===================================================================================================================================
!> Copies an existing state or timeavg file as a basis and write the merged data into it.
!> This ensures that all relevant information is contained in the merged file without too
!> much coding overhead.
!===================================================================================================================================
SUBROUTINE WriteTimeAverageByCopy(filename_in,filename_out,f,uavg)
! MODULES
USE MOD_HDF5_Output,ONLY:WriteArray,WriteAttribute
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN) :: filename_in  !< file to be copied
CHARACTER(LEN=*),INTENT(IN) :: filename_out !< output file
TYPE(tFileSet),INTENT(IN)   :: f            !< type with header and dataset information
REAL,INTENT(IN)             :: UAvg(f%totalsize) !< merged data (1D array with data for all datasets)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: offset2(maxDim) = 0
!===================================================================================================================================
CALL EXECUTE_COMMAND_LINE("cp -f "//TRIM(filename_in)//" "//TRIM(filename_out))

CALL OpenDataFile(filename_out,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
offset=0
DO i=1,f%nDataSets
  n=f%nDims(i)
  locsize(1:n)=f%nVal(1:n,i)
  startind=offset+1
  endind  =offset+PRODUCT(locsize(1:n))
  CALL WriteArray(TRIM(f%DataSetNames(i)),n,locsize(1:n),locsize(1:n),offset2(1:n),&
                          collective=.TRUE.,&
                          RealArray=UAvg(startInd:endInd))
  offset=endInd
END DO
CALL CloseDataFile()
END SUBROUTINE WriteTimeAverageByCopy


END PROGRAM TimeAvg
