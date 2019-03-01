#include "flexi.h"

!===================================================================================================================================
!> Add comments please!
!===================================================================================================================================
PROGRAM Nisp_RP
! MODULES
USE MOD_Globals
USE MOD_RPData                      ,ONLY:ReadRPData,AssembleRPData,FinalizeRPData
USE MOD_RPData_Vars                 
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_Nisp_RP_Vars
USE MOD_Nisp_RP_Output
USE MOD_IO_HDF5                     ,ONLY:DefineParametersIO_HDF5,InitIOHDF5
USE MOD_Spec                        ,ONLY:InitSpec,spec,FinalizeSpec
USE MOD_spec_Vars                        
USE MOD_RPSetVisuVisu_Vars                        
USE MOD_RPSetVisu                        
USE MOD_OutputRPVisu 
USE MOD_StringTools                 ,ONLY: STRICMP,GetFileExtension
USE MOD_RPInterpolation
USE MOD_RPInterpolation_Vars
USE MOD_ParametersVisu              ,ONLY: nVarVisu,VarNameVisu
USE MOD_RPSetVisuVisu_Vars          ,ONLY: nRP_global 
USE MOD_OutputRPVisu_Vars           ,ONLY:RPData_out,nSamples_out 
USE MOD_EquationRP
!#ifdef MPI
USE MOD_MPI                         ,ONLY:DefineParametersMPI,InitMPI
!#endif /* MPI */
USE MOD_Mesh                        ,ONLY: DefineParametersMesh,InitMesh
USE MOD_Mesh_Vars                   ,ONLY: nGlobalElems
USE MOD_Exactfunc                   ,ONLY: ExactFunc
USE MOD_Analyze                     ,ONLY: InitAnalyzeBasis
USE MOD_Output_Vars                 ,ONLY: ProjectName
USE MOD_EOS                         ,ONLY: ConsToPrim
USE MOD_EOS_Vars                    ,ONLY: KappaM1,R, Kappa
USE MOD_HDF5_Input                  ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataProps,ReadAttribute, ReadArray, GetDataSize
USE MOD_HDF5_Output                 ,ONLY: WriteAttribute,WriteArray,GenerateFileSkeleton,WriteHeader
USE MOD_Output_Vars                  ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output                  ,ONLY: insert_userblock
USE ISO_C_BINDING               ,ONLY: C_NULL_CHAR
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: iArg,iExt,nDataFiles,iSample, idx
REAL                               :: df
REAL                               :: dtSample = 1000, nSample_loc1,nSample_loc2
CHARACTER(LEN=255)                 :: InputIniFile
CHARACTER(LEN=255)                 :: InputDataFile
CHARACTER(LEN=255)                 :: sample_string,level_string,RPFilePath,iteration_string,iterationM1_string
CHARACTER(LEN=255),ALLOCATABLE     :: DataFiles(:)
INTEGER                              :: i, res, nVar_Field
CHARACTER(LEN=2)                     :: integer_string 
INTEGER                              :: line_no
LOGICAL                              :: FieldDataFound, create
LOGICAL                              :: userblockFound
CHARACTER(LEN=255)                   :: prmfile=".parameter.ini"
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Compute NISP modes based on RP data! ||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

CALL ParseCommandlineArguments()
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()


CALL DefineParameters()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
WRITE(UNIT_stdOut,'(A)') " || Compute NISP modes based on RP data!    ||"
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF

nArgs=COMMAND_ARGUMENT_COUNT()
IF(nArgs .LT. 2) CALL Abort(__STAMP__,'Missing argument')
CALL GET_COMMAND_ARGUMENT(1,InputIniFile)
! Get start index of file extension
iExt=INDEX(InputIniFile,'.',BACK = .TRUE.)

! check if first file is a .ini file
IF(InputIniFile(iExt+1:iExt+3) .NE. 'ini') &
CALL Abort(__STAMP__,'ERROR - No / invalid parameter file given.')
  
CALL prms%read_options(Args(1))



!======================================================
! ALLOCATE
!======================================================
nStochVar  = GETINT('NUncertainties')
N0         = GETINT('NStoch')
totSamples = GETINT('TotalSamples')
nDataFiles = nArgs-2
ALLOCATE(distributions(1:nStochVar))
ALLOCATE(weights(1:totSamples))
ALLOCATE(quadpoints(1:totSamples, 1:nStochVar))
ALLOCATE(DataFiles(1:nDataFiles))


!======================================================
! get list of input files
!======================================================
DO iArg=2,nArgs
  CALL GET_COMMAND_ARGUMENT(iArg,DataFiles(iArg-1))
END DO

! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: parameter_nisp.ini statefile_last_sample')
END IF


DO i=1,nStochVar
  write (integer_string,'(I0)') i
  distribution_string = disfunc//integer_string
  distributions(i)=  GETINT(distribution_string)
END DO

CALL InitParameters()
CALL InitIOHDF5()

!======================================================
! Read quadpoints and weights
!======================================================
OPEN(unit=99, FILE=file_weights,STATUS='old', ACTION='read')
DO i=1,totSamples
    READ(99,*)  weights(i)
END DO
CLOSE(99,STATUS='KEEP', IOSTAT=res)

OPEN(unit=99, FILE=file_quadpoints,STATUS='old', ACTION='read')
DO i=1,totSamples
    READ(99,*)  quadpoints(i,:)
END DO
CLOSE(99,STATUS='KEEP', IOSTAT=res)

!======================================================
! Initialize Multiindex and polynomials
!======================================================
CALL GetHermiteCoefficients()

CALL binom(nStochVar+N0,nStochVar, nPoly)

CALL CreateMultiIndex()

!======================================================
! timesample of samples
!======================================================
DO iSample=1,totSamples
   WRITE(sample_string,"(I0)") iSample
   ! readin RP Data from all input files
   DO iArg=1,nDataFiles
     InputDataFile=DataFiles(iArg)
     WRITE(UNIT_stdOut,'(132("="))')
     WRITE(UNIT_stdOut,'(A,I5,A,I5,A,I5)') ' PROCESSING FILE ',iArg,' of ',nDataFiles,' FILES'
     WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputDataFile),'" )'
     WRITE(UNIT_stdOut,'(132("="))')
     ! Get start index of file extension to check if it is a h5 file
     iExt=INDEX(InputDataFile,'.',BACK = .TRUE.)
     IF(InputDataFile(iExt+1:iExt+2) .NE. 'h5') &
       CALL Abort(__STAMP__,'ERROR - Invalid file extension!')
     ! Read in main attributes from given HDF5 State File
     InputDataFile='sample_'//TRIM(sample_string)//'/'//TRIM(InputDataFile)
     WRITE(UNIT_stdOUT,*) "READING DATA FROM RP FILE """,TRIM(InputDataFile), """"
     IF(iArg.EQ.1) THEN
       CALL ReadRPData(InputDataFile,firstFile=.TRUE.)
     ELSE 
       CALL ReadRPData(InputDataFile)
     END IF
   END DO
   
   CALL AssembleRPData()
   IF (iSample .EQ. 1) THEN
     CALL InitEquationRP()
     CALL InitInterpolation()
     CALL InitSpec()
     CALL InitOutput()
     
     ALLOCATE(UTimeseries(1:totSamples,1:nVarVisu,nRP_global,nSamples_out))
     ALLOCATE(UMeanTimeseries(1:nVarVisu,nRP_global,nSamples_out))
     ALLOCATE(UVarTimeseries(1:nVarVisu,nRP_global,nSamples_out))
     ALLOCATE(UModeTimeseries(1:nVarVisu,nRP_global,nSamples_out))
     UTimeseries =0.
     UMeanTimeseries =0.
     UVarTimeseries =0.
     UModeTimeseries = 0.
     ALLOCATE(time(nSamples_out))
   END IF

   CALL InterpolateEquiTime() !RPData
   CALL CalcEquationRP() !RPData_out
   uTimeseries(iSample,:,:,:)= RPData_out
   
   CALL spec()
   IF (iSample .EQ. 1) THEN
     ALLOCATE(UFFT(1:totSamples,1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UMeanFFT(1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UVarFFT(1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UModeFFT(1:nVarVisu,nRP_global,nSamples_spec))
     UFFT =0.
     UMeanFFT =0.
     UVarFFT =0.
     UModeFFT =0.
   END IF
     
   UFFT(iSample,:,:,:)=RPData_spec
   CALL WriteRPToHDF5(RPData_spec)
   IF(iSample.LT.totSamples) THEN
     CALL FinalizeRPSet()
     CALL FinalizeRPData()
   END IF
END DO
time = RPTime

!----------------------------------------------------------------------------------------------------------------------------------!
!EXPECTATION TIMESERIES & FFT
!----------------------------------------------------------------------------------------------------------------------------------!
CALL ComputeMode(0)
!----------------------------------------------------------------------------------------------------------------------------------!
!VARIANCE
!----------------------------------------------------------------------------------------------------------------------------------!
DO i=1, nPoly-1
  CALL ComputeMode(i)
  UVarTimeseries = UVarTimeseries + UModeTimeseries*UModeTimeseries
  UVarFFT = UVarFFT + UModeFFT*UModeFFT
END DO

SWRITE(UNIT_stdOut,'(A)') ' WRITING  MEAN and VARIANCE...'
CALL WriteMeanAndVarianceToHDF5()

CALL FinalizeInterpolation()
CALL FinalizeEquationRP()
CALL FinalizeSpec

SDEALLOCATE(UTimeseries)
SDEALLOCATE(UMeanTimeseries)
SDEALLOCATE(UVarTimeseries)
SDEALLOCATE(UModeTimeseries)
SDEALLOCATE(time)
SDEALLOCATE(UFFT)
SDEALLOCATE(UMeanFFT)
SDEALLOCATE(UVarFFT)
SDEALLOCATE(UModeFFT)
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif


END PROGRAM Nisp_RP

SUBROUTINE ComputeMode(i)
USE MOD_Nisp_RP_Vars
USE MOD_EOS                     ,ONLY: DefineParametersEOS,InitEOS, ConsToPrim
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataProps,ReadAttribute, ReadArray, GetDataSize
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER, INTENT(IN) :: i
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: j,k,l,m,iElem
REAL              :: evalPoly, y_out, y_out_dummy
!-----------------------------------------------------------------------------------------------------------------------------------
IF(i==0) THEN
  DO j=1, totSamples
    UMeanTimeseries = UMeanTimeseries + uTimeseries(j,:,:,:)*weights(j)
    UMeanFFT        = UMeanFFT        + UFFT(j,:,:,:)*weights(j)
  END DO
ELSE
  y_out=0.
  UModeTimeseries = 0.
  UModeFFT = 0.
  DO j=1, totSamples
    evalPoly = 1.
    DO k=1, nStochVar
      IF(distributions(k)==1) THEN !Legendre on [-0.5,0.5]
        CALL LegendrePolynomialAndDerivative(MultiIndex(i,k), 2*quadpoints(j,k)-1.,.FALSE.,y_out,y_out_dummy)
      ELSE !Hermite
        CALL EvaluateHermitePoly(MultiIndex(i,k),quadpoints(j,k),y_out)
      END IF
      evalPoly = evalPoly*y_out
    END DO
    UModeTimeseries = UModeTimeseries+ uTimeseries(j,:,:,:)*evalPoly*weights(j)
    UModeFFT = UModeFFT+ UFFT(j,:,:,:)*evalPoly*weights(j)
  END DO  
END IF


END SUBROUTINE
!===================================================================================================================================
!> Initialize parameter variables of Posti tool
!===================================================================================================================================
SUBROUTINE DefineParameters()
! MODULES
USE MOD_ReadInTools 
USE MOD_Nisp_RP_Vars,             ONLY: disfunc

IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)                   :: distribution_string
INTEGER                              :: i
CHARACTER(LEN=3)                     :: integer_string 
!===================================================================================================================================

CALL prms%SetSection('Nisp Parameters for RP')
CALL prms%CreateStringOption(   "ProjectName"        , "Define Name of the project.",multiple=.TRUE.)
CALL prms%CreateStringOption(   "RP_DefFile"         , "File which defines the RP setup.",multiple=.TRUE.)


CALL prms%CreateLogicalOption('OutputPoints'       ,"TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "doPSD"              , "Perform discrete fourier transform",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "UseNonDimensionalEqn"              , "Perform discrete fourier transform",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "SkipSample"         , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "nBlocks"            , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "BlockSize"          , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "RP_specified"       , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "doFFT"              , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "hanning"            , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  'fourthDeriv'        , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  'ThirdOct'           , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'SamplingFreq'       , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'cutoffFreq'         , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'u_inf'              , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'mu0'                , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'chord'              , "TODO",multiple=.TRUE.)
CALL prms%CreateStringOption(   'VarName'            ,"TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    ('OutputFormat'       ,"TODO",multiple=.TRUE.)
CALL prms%SetSection("NISP Setup")
CALL prms%CreateIntOption(      "TotalSamples"             , "Number of samples which need to be computedl",multiple=.TRUE.)

CALL prms%SetSection("Stochastic Input")
CALL prms%CreateIntOption( "NUncertainties"             , "Number of stochastic input variables", multiple=.TRUE.)
CALL prms%CreateIntOption( "NStoch"             , "max. polynomial chaos degree",multiple=.TRUE.)
CALL prms%SetSection("Equation of State")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent",multiple=.TRUE.)
CALL prms%CreateRealOption(     'R',            "Specific gas constant",multiple=.TRUE.)
CALL prms%CreateRealOption(     'Pr',            "Prandtl number",multiple=.TRUE.)
CALL prms%CreateRealOption(     'mu0',            "Dynamic Viscosity",multiple=.TRUE.)
DO i=1,100
  write (integer_string,'(I0)') i
  distribution_string = disfunc//integer_string
  CALL prms%CreateIntOption(      distribution_string             , "1: uniform distribution  2: Normal distributuion")
END DO
END SUBROUTINE DefineParameters


!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Globals
USE MOD_Readintools         
USE MOD_Nisp_RP_Vars
USE MOD_ParametersVisu
USE MOD_EOS_Vars                ,ONLY: KappaM1,R, Kappa
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

SUBROUTINE CreateMultiIndex()
! INPUT / OUTPUT VARIABLES 
!N0 = highest polynomial degree |i|=N0, i=(i1,...,i_stochDim)
!nPoly = #number of stochdim-dimensional chaos polynomials of order <=N0
!Algorithm in LeMaitre, p.524
!MODULES
USE MOD_Nisp_RP_Vars, ONLY: MultiIndex, nStochVar, N0, nPoly
IMPLICIT NONE!
!----------------------------------------------------------------------------------------------------------------------------------!

! LOCAL VARIABLES
INTEGER           :: i, j, k, m, P, L
INTEGER           :: p2(1:N0,1:nStochVar)
!----------------------------------------------------------------------------------------------------------------------------------!

ALLOCATE(MultiIndex(0:nPoly,1:nStochVar))
MultiIndex = 0.
DO i=1,nStochVar
  MultiIndex(i,i) = 1
END DO
p2= 0.
P= nStochVar
p2(1,:)=1
DO k=2,N0
  L=P
  DO i=1,nStochVar
    DO m=i,nStochVar
      p2(k,i)= p2(k,i) + p2(k-1,m)
    END DO
  END DO
  DO j=1,nStochVar
    DO m= L- p2(k,j)+1, L
      P=P+1
      MultiIndex(P,:)=MultiIndex(m,:)
      MultiIndex(P,j)=MultiIndex(P,j)+1
    END DO
  END DO
END DO

END SUBROUTINE CreateMultiIndex

!===================================================================================================================================
!> Evaluate the Legendre polynomial L_N and its derivative at position x[-1,1]
!> recursive algorithm using the N_in-1 N_in-2 Legendre polynomials
!> algorithm 22, Kopriva book
!===================================================================================================================================
SUBROUTINE LegendrePolynomialAndDerivative(N_in,xIn,HalfInterval,L,Lder) !YYY
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in         !< (IN)  polynomial degree, (N+1) CLpoints
DOUBLE PRECISION,INTENT(IN)    :: xIn          !< (IN)  coordinate value in the interval [-1,1] !YYY
LOGICAL,INTENT(IN) :: HalfInterval !< (IN)  input coordinates are for interval [-0.5,0.5] !YYY
DOUBLE PRECISION,INTENT(OUT)   :: L            !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
DOUBLE PRECISION,INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iLegendre
DOUBLE PRECISION    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
DOUBLE PRECISION    :: x !YYY
DOUBLE PRECISION    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF (HalfInterval) THEN !YYY
  x=2.*xIn !YYY
ELSE !YYY
  x=xIn !YYY
END IF !YYY
IF(N_in .EQ. 0)THEN
  L=1.
  Lder=0.
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.
ELSE ! N_in > 1
  L_Nm2=1.
  L_Nm1=x
  Lder_Nm2=0.
  Lder_Nm1=1.
  DO iLegendre=2,N_in
    L=(REAL(2*iLegendre-1)*x*L_Nm1 - REAL(iLegendre-1)*L_Nm2)/REAL(iLegendre)
    Lder=Lder_Nm2 + REAL(2*iLegendre-1)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize Polynomials for new Interval !YYY
IF (HalfInterval) THEN !YYY
  L=L*SQRT(REAL(2*N_in)+1.) !YYY
  Lder=Lder*SQRT(REAL(2*N_in)+1.) !YYY
ELSE !YYY
  L=L*SQRT(REAL(N_in)+0.5)
  Lder=Lder*SQRT(REAL(N_in)+0.5)
END IF !YYY
END SUBROUTINE LegendrePolynomialAndDerivative




SUBROUTINE GetHermiteCoefficients()
! MODULES
USE MOD_Nisp_RP_Vars, ONLY:HermiteCoeff,N0
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p
!==================================================================================================================================
ALLOCATE(HermiteCoeff(0:N0,0:N0))

HermiteCoeff(:,:)=0.
HermiteCoeff(0,0)=1.
HermiteCoeff(1,1)=1.
DO p=2,N0
    HermiteCoeff(1:p,p)=HermiteCoeff(0:p-1,p-1)/SQRT(REAL(p))
    HermiteCoeff(0:p-2,p)=HermiteCoeff(0:p-2,p)-HermiteCoeff(0:p-2,p-2)*SQRT(REAL(p-1)/REAL(p))
END DO
END SUBROUTINE GetHermiteCoefficients


!===================================================================================================================================
!> evaluate Hermite polynomial of degree pStoch at xEval
!===================================================================================================================================
SUBROUTINE EvaluateHermitePoly(pStoch,xEval,y_out)
! MODULES
USE MOD_Nisp_RP_Vars, ONLY:HermiteCoeff 
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: pStoch  !< (IN)  input polynomial degree
REAL, INTENT(IN)    :: xEval
REAL, INTENT(OUT)   :: y_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p
!==================================================================================================================================
y_out=0.
DO p=0,pStoch
  y_out=y_out+HermiteCoeff(p,pStoch)*xEval**p  
END DO
END SUBROUTINE EvaluateHermitePoly


SUBROUTINE binom(n, k, resu)
! MODULES
IMPLICIT NONE!
! INPUT / OUTPUT VARIABLES 
INTEGER , INTENT(IN)   :: n,k
INTEGER , INTENT(OUT)  :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j
j=k
resu=1
IF(k.GT.(n-k)) THEN
  j = n - k;
END IF
DO i=0,j-1
  resu = resu*(n - i)
  resu = resu/(i + 1);
END DO
END SUBROUTINE binom

