#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM EstimateSigma_RP
! MODULES
USE MOD_Globals
USE MOD_RPData                      ,ONLY:ReadRPData,AssembleRPData,FinalizeRPData
USE MOD_RPData_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools                 ,ONLY:PrintDefaultParameterFile
USE MOD_ReadInTools
USE MOD_EstimateSigma_RP_Vars
USE MOD_EstimateSigma_RP_Output
USE MOD_EstimateSigma_RP_Input
USE MOD_IO_HDF5                     ,ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Spec                        ,ONLY:InitSpec,spec,FinalizeSpec
USE MOD_spec_Vars
USE MOD_RPSetVisuVisu_Vars
USE MOD_RPSetVisu
USE MOD_OutputRPVisu
USE MOD_RPInterpolation
USE MOD_RPInterpolation_Vars        ,ONLY:CalcTimeAverage,dt
USE MOD_ParametersVisu              ,ONLY: nVarVisu,VarNameVisu
USE MOD_EquationRP
USE MOD_EOS                         ,ONLY:DefineParametersEOS,InitEOS
!#ifdef MPI
USE MOD_MPI                         ,ONLY:DefineParametersMPI,InitMPI
!#endif /* MPI */
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iExt,nDataFiles,iSample
REAL                               :: df
CHARACTER(LEN=255)                 :: InputIniFile
CHARACTER(LEN=255)                 :: InputDataFile
CHARACTER(LEN=255)                 :: sample_string,level_string,RPFilePath,iteration_string,iterationM1_string
CHARACTER(LEN=255),ALLOCATABLE     :: DataFiles(:)
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Compute SigmaSq based on RP data!    ||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

CALL ParseCommandlineArguments()
CALL DefineParameters()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersEOS()


! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
WRITE(UNIT_stdOut,'(A)') " || Compute sigmaSq based on RP data!    ||"
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs .LT. 2) CALL Abort(__STAMP__,'Missirg argument')
CALL GET_COMMAND_ARGUMENT(1,InputIniFile)
! Get start index of file extension
iExt=INDEX(InputIniFile,'.',BACK = .TRUE.)
! check if first file is a .ini file
IF(InputIniFile(iExt+1:iExt+3) .NE. 'ini') &
  CALL Abort(__STAMP__,'ERROR - No / invalid parameter file given.')

CALL prms%read_options(Args(1))

!CALL InitRPSet()

nDataFiles=nArgs-1
ALLOCATE(DataFiles(1:nDataFiles))

! get list of input files
DO iArg=2,nArgs
  CALL GET_COMMAND_ARGUMENT(iArg,DataFiles(iArg-1))
END DO

CALL InitParameters()
CALL InitIOHDF5()
!-----------------------------------------------------------------------------------------------------------------------------------
WRITE(level_string,"(I0)") iLevel
WRITE(iteration_string,"(I0)") nIter
WRITE(iterationM1_string,"(I0)") nIterIn
Filename_SigmaSq  = 'level_'//TRIM(level_string)//'/'//'sigmasq_spec.dat'
Filename_SigmaSq_Fine  = 'level_'//TRIM(level_string)//'/'//'fine_sigmasq_spec.dat'
FileNameSumsIn    = 'level_'//TRIM(level_string)//'/sums_'//TRIM(iterationM1_string)//'_spec.h5'
FileNameSumsOut   = 'level_'//TRIM(level_string)//'/sums_'//TRIM(iteration_string)//'_spec.h5'
FilenameMean      = 'level_'//TRIM(level_string)//'/'//'mean_spec.h5'
FilenameVariance  = 'level_'//TRIM(level_string)//'/'//'variance_spec.h5'
snSamples  =1./REAL(nEnd)
snSamplesM1=1./(REAL(nEnd)-1.)

DO iSample=nStart,nEnd
   WRITE(sample_string,"(I0)") iSample

   ! readin RP Data from all input files
   DO iArg=1,nDataFiles
     InputDataFile=DataFiles(iArg)
     WRITE(UNIT_stdOut,'(132("="))')
     WRITE(UNIT_stdOut,'(A,I5,A,I5,A,I5)') ' PROCESSING FILE ',iArg,' of ',nDataFiles,' FILES on the fine grid of Level.',iLevel
     WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputDataFile),'" )'
     WRITE(UNIT_stdOut,'(132("="))')


     ! Get start index of file extension to check if it is a h5 file
     iExt=INDEX(InputDataFile,'.',BACK = .TRUE.)
     IF(InputDataFile(iExt+1:iExt+2) .NE. 'h5') &
       CALL Abort(__STAMP__,'ERROR - Invalid file extension!')
     ! Read in main attributes from given HDF5 State File
     InputDataFile='level_'//TRIM(level_string)//'/f/sample_'//TRIM(sample_string)//'/'//&
      TRIM(InputDataFile)
     WRITE(UNIT_stdOUT,*) "READING DATA FROM RP FILE """,TRIM(InputDataFile), """"
     IF(iArg.EQ.1) THEN
       CALL ReadRPData(InputDataFile,firstFile=.TRUE.)
     ELSE
       CALL ReadRPData(InputDataFile)
     END IF
   END DO

   ! assemble data to one global array
   CALL AssembleRPData()
   IF (iSample .EQ. nStart) THEN
     CALL InitEquationRP()
     CALL InitInterpolation()
     CALL InitSpec()
     CALL InitOutput()
   END IF
   !CALL InitOutput()
   CALL InterpolateEquiTime()
   CALL CalcEquationRP()
   CALL spec()


   IF (iSample .EQ. nStart) THEN
     !ALLOCATE(UFine(nVarVisu,nPoints,nSamples_spec))
     !ALLOCATE(UCoarse(nVarVisu,nPoints,nSamples_spec))
     !ALLOCATE(DeltaU(nVarVisu,nPoints,nSamples_spec))
     !ALLOCATE(USum(nVarVisu,nPoints,nSamples_spec))
     !ALLOCATE(USqSum(nVarVisu,nPoints,nSamples_spec))
     !ALLOCATE(SigmaSqSpec(nVarVisu,nPoints,nSamples_spec))
     !UFine(:,:,:)=0
     !UCoarse(:,:,:)=0
     !DeltaU(:,:,:)=0
     !USum(:,:,:)=0
     !USqSum(:,:,:)=0
     ALLOCATE(UFine(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(UCoarse(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(UFineSum(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(UCoarseSum(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(UFineSqSum(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(UCoarseSqSum(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(DUSqSum(nVarVisu,nPoints,nSamples_spec))
     ALLOCATE(SigmaSqSpec(nVarVisu,nPoints,nSamples_spec))
     UFine(:,:,:)=0
     UCoarse(:,:,:)=0
     UFineSum(:,:,:)=0
     UCoarseSum(:,:,:)=0
     UFineSqSum(:,:,:)=0
     UCoarseSqSum(:,:,:)=0
     DUSqSum(:,:,:)=0
     SigmaSqSpec(:,:,:)=0
     IF (.NOT. (nStart.EQ.1)) THEN
       CALL ReadSumsFromHDF5()
     END IF
   END IF

   UFine(:,:,:)=RPData_spec(:,:,:)
   CALL FinalizeRPData()
   !CALL FinalizeOutput()
   IF (iLevel .GT. 1 .OR. iSample .LT. nEnd) CALL FinalizeSpec()
   IF (iLevel .GT. 1 .OR. iSample .LT. nEnd) CALL FinalizeRPSet()


   IF (iLevel .GT. 1) THEN
     ! readin RP Data from all input files
     DO iArg=1,nDataFiles
       InputDataFile=DataFiles(iArg)
       WRITE(UNIT_stdOut,'(132("="))')
       WRITE(UNIT_stdOut,'(A,I5,A,I5,A,I5)') ' PROCESSING FILE ',iArg,' of ',nDataFiles,' FILES on the coarse grid of Level.',iLevel
       WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputDataFile),'" )'
       WRITE(UNIT_stdOut,'(132("="))')


       ! Get start index of file extension to check if it is a h5 file
       iExt=INDEX(InputDataFile,'.',BACK = .TRUE.)
       IF(InputDataFile(iExt+1:iExt+2) .NE. 'h5') &
         CALL Abort(__STAMP__,'ERROR - Invalid file extension!')
       ! Read in main attributes from given HDF5 State File
       InputDataFile='level_'//TRIM(level_string)//'/c/sample_'//TRIM(sample_string)//'/'//&
        TRIM(InputDataFile)
       WRITE(UNIT_stdOUT,*) "READING DATA FROM RP FILE """,TRIM(InputDataFile), """"
       IF(iArg.EQ.1) THEN
         CALL ReadRPData(InputDataFile,firstFile=.TRUE.)
       ELSE
         CALL ReadRPData(InputDataFile)
       END IF
     END DO

     ! assemble data to one global array
     CALL AssembleRPData()
     !CALL InitOutput()
     CALL InterpolateEquiTime()
     CALL CalcEquationRP()
     CALL spec()

     UCoarse(:,:,:)=RPData_spec(:,:,:)

     CALL FinalizeRPData()
     !CALL FinalizeOutput()
     IF (iSample .LT. nEnd) CALL FinalizeSpec()
     IF (iSample .LT. nEnd) CALL FinalizeRPSet()
   END IF


   UFineSum     = UFineSum     + UFine
   UCoarseSum   = UCoarseSum   + UCoarse

   UFineSqSum   = UFineSqSum   + UFine*UFine
   UCoarseSqSum = UCoarseSqSum + UCoarse*UCoarse

   DUSqSum      = DUSqSum      + (UFine-UCoarse)*(UFine-UCoarse)


END DO

! SigmaSq Computation
df=RPData_freq(nSamples_spec)/(nSamples_Spec-1)



CALL WriteSumsToHDF5()

SigmaSqSpec = snSamplesM1*( DUSqSum - snSamples*(UFineSum-UCoarseSum) * (UFineSum - UCoarseSum) )
SigmaSq = Sum(df*SigmaSqSpec(VarAna,RP_specified,3:nSamples_Spec))
! Output SigmaSq
SWRITE(UNIT_stdOut,'(132("="))')
CALL WriteSigmaSq(Filename_SigmaSq)

! Write SigmaSq only on fine grid
SigmaSqSpec = snSamplesM1*( UFineSqSum - snSamples * UFineSum * UFineSum )
SigmaSq = Sum(df*SigmaSqSpec(VarAna,RP_specified,3:nSamples_Spec))
SWRITE(UNIT_stdOut,'(132("="))')
CALL WriteSigmaSq(Filename_SigmaSq_Fine)

SigmaSqSpec = SigmaSqSpec - snSamplesM1*( UCoarseSqSum - snSamples * UCoarseSum*UCoarseSum )
CALL WriteMeanAndVarianceToHDF5()

CALL FinalizeRPSet()
CALL FinalizeSpec()
CALL FinalizeOutput()
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A,I2,A)') 'SigmaSq on Level ', iLevel,' computed !'
WRITE(UNIT_stdOut,'(A,ES14.7)') 'SigmaSq: ', SigmaSq
WRITE(UNIT_stdOut,'(A,I1,A)') 'SigmaSq written to level_', iLevel,'/SigmaSq.dat !'
WRITE(UNIT_stdOut,'(132("="))')

CALL FinalizeParameters()
CALL FinalizeEquationRP()
CALL FinalizeInterpolation()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') 'SigmaSq comptued from RECORDPOINTS done ! '
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

CALL prms%CreateIntOption(      "nStart"             , "Number of first new sample")
CALL prms%CreateIntOption(      "nEnd"               , "Number of last new sample")
CALL prms%CreateIntOption(      "nIter"              , "Number of current iteration")
CALL prms%CreateIntOption(      "nIterIn"              , "Number of current iteration")
CALL prms%CreateIntOption(      "iLevel"             , "SigmaSq is computed only on this specific level")
CALL prms%CreateIntOption(      "varAna"             , "Specifiy the variable from U cons to comput SigmaSq")

CALL prms%CreateLogicalOption('OutputPoints'       ,"TODO",".FALSE.")
CALL prms%CreateLogicalOption(  "doPSD"              , "Perform discrete fourier transform",".TRUE.")
CALL prms%CreateLogicalOption(  "UseNonDimensionalEqn"              , "Perform discrete fourier transform",".TRUE.")
CALL prms%CreateIntOption    (  "SkipSample"         , "TODO")
CALL prms%CreateIntOption    (  "nBlocks"            , "TODO")
CALL prms%CreateIntOption    (  "BlockSize"          , "TODO")
CALL prms%CreateIntOption    (  "RP_specified"       , "TODO")
CALL prms%CreateLogicalOption(  "doFFT"              , "TODO",".FALSE.")
CALL prms%CreateLogicalOption(  "hanning"            , "TODO",".FALSE.")
CALL prms%CreateLogicalOption(  'fourthDeriv'        , "TODO",".FALSE.")
CALL prms%CreateLogicalOption(  'ThirdOct'           , "TODO",".FALSE.")
CALL prms%CreateRealOption   (  'SamplingFreq'       , "TODO")
CALL prms%CreateRealOption   (  'cutoffFreq'         , "TODO")
CALL prms%CreateRealOption   (  'u_inf'              , "TODO")
CALL prms%CreateRealOption   (  'mu0'                , "TODO")
CALL prms%CreateRealOption   (  'R'                  , "TODO")
CALL prms%CreateRealOption   (  'kappa'              , "TODO")
CALL prms%CreateRealOption   (  'Pr'                 , "TODO")
CALL prms%CreateRealOption   (  'chord'              , "TODO")
CALL prms%CreateStringOption(   'VarName'            ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    ('OutputFormat'       ,"TODO")

END SUBROUTINE DefineParameters


!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Globals
USE MOD_Readintools         ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_EstimateSigma_RP_Vars
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
iLevel=GETINT('iLevel')
nStart=GETINT('nStart')
nEnd=GETINT('nEnd')
nIter=GETINT('nIter')
nIter=GETINT('nIterIn')
VarAna=GETINT('VarAna')
RP_specified=GETINT('RP_specified')

OutputPoints    = GETLOGICAL('OutputPoints','.TRUE.')
OutputFormat    = GETINT('OutputFormat','2')


! =============================================================================== !
! FOURIER TRANSFORM
! =============================================================================== !
skip = GETINT('SkipSample','1')
doFFT=GETLOGICAL('doFFT','F')
doPSD=GETLOGICAL('doPSD','T')
IF(doFFT.OR.doPSD) doSpec=.TRUE.
IF(doSpec) THEN
  ! two readin "modes" for spectrum averaging:
  ! 1. Prescription of number of blocks
  nBlocks=GETINT('nBlocks','1')
  ! 2. Prescription of Sampling Frequency and Blocksize
  samplingFreq=GETREAL('SamplingFreq','-999')
  IF(samplingFreq.GT.0.) THEN
    BlockSize=GETINT('BlockSize')
  END IF
  doHanning=GETLOGICAL('hanning','F')
  fourthDeriv=GETLOGICAL('fourthDeriv','F')
  ThirdOct=GETLOGICAL('ThirdOct','F')
  IF (ThirdOct) THEN
    u_infPhys   = GETREAL('u_inf') !velocity for re-dimensionalization of frequency
    chordPhys = GETREAL('chord')   !length for re-dimensionalization of frequency
  END IF
END IF

IF(doSpec) cutoffFreq=GETREAL('cutoffFreq','-999.')

equiTimeSpacing=.TRUE.
END SUBROUTINE InitParameters


!===================================================================================================================================
!> This routine builds the mappings from the total number of variables available for visualization to number of calculation
!> and visualization variables.
!>  1. Read 'VarName' options from the parameter file. This are the quantities that will be visualized.
!>  2. Initialize the dependecy table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapVisu' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!===================================================================================================================================
SUBROUTINE Build_mapCalc_mapVisu()
USE MOD_ParametersVisu
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
ALLOCATE(mapVisu(1:nVarDep))
mapVisu = 0
nVarVisuTotal = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,nVarVisu
  DO iVar2=1,nVarDep
    IF (STRICMP(VarNameVisu(iVar), VarNamesAll(iVar2))) THEN
      mapVisu(iVar2) = nVarVisuTotal+1
      nVarVisuTotal = nVarVisuTotal + 1
    END IF
  END DO
END DO

! Calculate all dependencies:
! For each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarDep
  DepTable(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable(iVar,iVar2).EQ.1) &
      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
  END DO
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
ALLOCATE(mapCalc(1:nVarDep))
mapCalc = 0
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    mapCalc = MAX(mapCalc,DepTable(iVar,1:nVarDep))
  END IF
END DO
! enumerate mapCalc
nVarCalc = 0
DO iVar=1,nVarDep
  IF (mapCalc(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc(iVar) = nVarCalc
  END IF
END DO

! print the dependecy table
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarDep
  WRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarNamesAll(iVar))
END DO

! print the mappings
WRITE(format,'(I2)') nVarDep
WRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
WRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

END SUBROUTINE Build_mapCalc_mapVisu
