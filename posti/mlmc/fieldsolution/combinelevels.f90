#include "flexi.h"

!===================================================================================================================================
!> This tool will take pre-averaged files (TimeAvg, Flucs) or simple state files 
!> and perform global temporal averaging
!===================================================================================================================================
PROGRAM Combinelevels
! MODULES
USE MOD_Globals
USE MOD_ReadInTools, ONLY:prms
USE MOD_Commandline_Arguments
USE MOD_StringTools,ONLY:STRICMP
USE MOD_IO_HDF5
USE MOD_MPI,        ONLY:InitMPI
USE MOD_MLMC_Vars
USE MOD_MLMC_Input, ONLY:GetParams,ReadSums
USE MOD_MLMC_Output, ONLY:WriteMeanAndVarianceToHDF5
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! TYPE AND PARAMETER DEFINITIONS
INTEGER             :: iArg
INTEGER             :: stat,iniUnit
CHARACTER(LEN=255)  :: FileNameDummy=".dummy.ini"
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
OPEN(NEWUNIT= iniUnit,        &
     FILE   = FileNameDummy, &
     STATUS = 'NEW',          &
     IOSTAT = stat)
CLOSE(iniUnit)
CALL prms%read_options(FileNameDummy)
CALL DefineParametersIO_HDF5()
CALL InitIOHDF5()
CALL ParseCommandlineArguments()
IF ((doPrintHelp.GT.0) .OR. (nArgs.LT.1)) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: ./mlmc_combinelevels [sums_xy.h5 ...]')
END IF

CALL GetParams(TRIM(Args(1)))
ALLOCATE(UFineSum    (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(UCoarseSum  (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(UFineSqSum  (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(UCoarseSqSum(nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(DUSqSum     (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(Mean        (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
ALLOCATE(StdDev      (nVal(1),nVal(2),nVal(3),nVal(4),nVal(5)))
Mean=0.
StdDev=0.

! Start the averaging
WRITE(UNIT_stdOut,'(132("="))')

DO iArg=1,nArgs
  SWRITE(UNIT_stdOut,'(A,A,A)') "Processing file ",Args(iArg),"..."
  CALL ReadSums(Args(iArg))
  snSamples=1./nSamples_Sums
  snSamplesM1=1./(nSamples_Sums-1)

  Mean = Mean + snSamples*(UFineSum-UCoarseSum)
  StdDev = StdDev + snSamplesM1*((UFineSqSum-UCoarseSum)   - snSamples*(UFineSum**2-UCoarseSum**2))
END DO
StdDev = SQRT(ABS(StdDev))

CALL WriteMeanAndVarianceToHDF5()

SDEALLOCATE(UFineSum)
SDEALLOCATE(UCoarseSum)
SDEALLOCATE(UFineSqSum)
SDEALLOCATE(UCoarseSqSum)
SDEALLOCATE(DUSqSum)
SDEALLOCATE(Mean)
SDEALLOCATE(StdDev)
 
SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') "Merging DONE!"
SWRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Combinelevels
