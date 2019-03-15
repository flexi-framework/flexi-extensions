#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM EstimateSigma_RP
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
!USE MOD_ReadInTools                 ,ONLY:PrintDefaultParameterFile
USE MOD_ReadInTools
USE MOD_CombineLevels_RP_Vars
USE MOD_IO_HDF5                     ,ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_IO_HDF5                     ,ONLY:InitMPIInfo
!USE MOD_RPSetVisuVisu_Vars                        
USE MOD_RPSetVisu                        
!USE MOD_OutputRPVisu                      
USE MOD_CombineLevels_RP_Output
USE MOD_CombineLevels_RP_Input
USE MOD_ParametersVisu              ,ONLY: RP_DefFile
!#ifdef MPI
USE MOD_MPI                         ,ONLY:DefineParametersMPI,InitMPI
!#endif /* MPI */
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iExt,nDataFiles,iLevel
REAL                               :: df
CHARACTER(LEN=255)                 :: InputIniFile
CHARACTER(LEN=255)                 :: level_string,Iter_string,RPFilePath,iteration_string,iterationM1_string
CHARACTER(LEN=255),ALLOCATABLE     :: DataFiles(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
CALL InitMPIInfo()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Combine Levels based on RP data!    ||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

CALL ParseCommandlineArguments()
CALL DefineParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
!CALL DefineParametersEOS()


! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
WRITE(UNIT_stdOut,'(A)') " || Compute sigmaSq based on RP data!    ||"
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs .LT. 1) CALL Abort(__STAMP__,'Missirg argument')
CALL GET_COMMAND_ARGUMENT(1,InputIniFile)
! Get start index of file extension
iExt=INDEX(InputIniFile,'.',BACK = .TRUE.)
! check if first file is a .ini file
IF(InputIniFile(iExt+1:iExt+3) .NE. 'ini') &
  CALL Abort(__STAMP__,'ERROR - No / invalid parameter file given.')

CALL prms%read_options(Args(1))

CALL InitParameters()
CALL InitIOHDF5()

CALL InitRPSet(RP_DefFile)
!-----------------------------------------------------------------------------------------------------------------------------------

DO iLevel=1,nLevels
   WRITE(level_string,"(I0)") iLevel

   WRITE(UNIT_stdOut,'(132("="))')
   WRITE(UNIT_stdOut,'(A,I5,A,I5,A,I5)') ' PROCESSING Level ',iLevel,' of ',nLevels,'Levels.'
   WRITE(UNIT_stdOut,'(132("="))')
  

   ! Get start index of file extension to check if it is a h5 file
   !WRITE(iter_string,"(I0)") nIter
   FileNameMean='level_'//TRIM(level_string)//'/mean_spec.h5'
   iExt=INDEX(FileNameMean,'.',BACK = .TRUE.)
   IF(FileNameMean(iExt+1:iExt+2) .NE. 'h5') &
     CALL Abort(__STAMP__,'ERROR - Invalid file extension!')


   FileNameVariance='level_'//TRIM(level_string)//'/variance_spec.h5'
   iExt=INDEX(FileNameVariance,'.',BACK = .TRUE.)
   IF(FileNameVariance(iExt+1:iExt+2) .NE. 'h5') &
     CALL Abort(__STAMP__,'ERROR - Invalid file extension!')


   ! Read in main attributes from given HDF5 State File
   !WRITE(UNIT_stdOUT,*) "READING DATA FROM RP FILE """,TRIM(FileNameMean), """"

   IF (iLevel .EQ. 1) CALL InitReadSumsFromHDF5()
   CALL ReadSumsFromHDF5()

END DO



SWRITE(UNIT_stdOut,'(132("="))')
CALL WriteMeanAndVarianceToHDF5()

CALL FinalizeReadSumsFromHDF5()
CALL FinalizeRPSet()
!CALL FinalizeOutput()

CALL FinalizeParameters()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') 'Final Mean and Variance comptued from RECORDPOINTS done ! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM EstimateSigma_RP


!===================================================================================================================================
!> Initialize parameter variables of Posti tool
!===================================================================================================================================
SUBROUTINE DefineParameters()
! MODULES
USE MOD_ReadInTools 
IMPLICIT NONE
!===================================================================================================================================

CALL prms%SetSection('MLMC Parameters for RP')
CALL prms%CreateStringOption(   "ProjectName"        , "Define Name of the project.")
CALL prms%CreateStringOption(   "RP_DefFile"         , "File which defines the RP setup.")
CALL prms%CreateIntOption(      "nLevels"               , "Number of last new sample")
!CALL prms%CreateIntOption(      "nIter"              , "Number of current iteration")
CALL prms%CreateIntOption(      "OutputFormat"       , "TODO")
CALL prms%CreateLogicalOption(      "OutputPoints"       , "TODO")

END SUBROUTINE DefineParameters


!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Globals
USE MOD_Readintools         ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_Combinelevels_RP_Vars
USE MOD_ParametersVisu
!USE MOD_RPInterpolation_Vars,ONLY:calcTimeAverage
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================

Projectname=GETSTR('ProjectName','')
RP_DefFile=GETSTR('RP_DefFile','')
! =============================================================================== !
! RP INFO
! =============================================================================== !
nLevels=GETINT('nLevels')
!nIter=GETINT('nIter')

OutputPoints    = GETLOGICAL('OutputPoints','.TRUE.')
OutputFormat    = GETINT('OutputFormat','2')

END SUBROUTINE InitParameters
