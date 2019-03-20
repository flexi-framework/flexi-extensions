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

MODULE MOD_NISP_RP

INTERFACE computeModes
  MODULE PROCEDURE computeModes
END INTERFACE

INTERFACE FinalizeNisp_RP
  MODULE PROCEDURE FinalizeNisp_RP
END INTERFACE

PUBLIC::DefineParametersNisp_RP,InitNisp_RP,PerformSampleFFT,computeModes,FinalizeNisp_RP
CONTAINS

!===================================================================================================================================
!> Define parameters of Posti Nisp RP tool
!===================================================================================================================================
SUBROUTINE DefineParametersNisp_RP()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ReadInTools 
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%SetSection('Nisp Parameters for RP')
CALL prms%CreateStringOption(   "RP_DefFile"         , "File which defines the RP setup.",multiple=.TRUE.)

CALL prms%CreateLogicalOption(  'OutputPoints'       , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "doPSD"              , "Perform discrete fourier transform",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "UseNonDimensionalEqn"              , "Perform discrete fourier transform",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "SkipSample"         , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "nBlocks"            , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  "BlockSize"          , "TODO",multiple=.TRUE.)
!CALL prms%CreateIntOption    (  "RP_specified"       , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "doFFT"              , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  "hanning"            , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  'fourthDeriv'        , "TODO",multiple=.TRUE.)
CALL prms%CreateLogicalOption(  'ThirdOct'           , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'SamplingFreq'       , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'cutoffFreq'         , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'u_inf'              , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'chord'              , "TODO",multiple=.TRUE.)
CALL prms%CreateStringOption (  'VarName'            , "TODO",multiple=.TRUE.)
CALL prms%CreateIntOption    (  'OutputFormat'       , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'kappa'              , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'Pr'                 , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'R'                  , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'mu0'                , "TODO",multiple=.TRUE.)
END SUBROUTINE DefineParametersNisp_RP

!===================================================================================================================================
!> Read parameters of Posti Nisp RP tool
!===================================================================================================================================
SUBROUTINE InitNisp_RP()
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Readintools         
USE MOD_Nisp_RP_Vars
USE MOD_ParametersVisu
USE MOD_IO_HDF5           ,ONLY: File_ID,OpenDataFile,CloseDataFile
USE MOD_HDF5_Input        ,ONLY: ReadAttribute,ReadArray
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)      :: StochFile
INTEGER                 :: iArg
!===================================================================================================================================
! =============================================================================== !
! Read in parameter.ini
! =============================================================================== !
RP_DefFile   = GETSTR('RP_DefFile','')
!RP_specified = GETINT('RP_specified')
OutputPoints = GETLOGICAL('OutputPoints','.TRUE.')
OutputFormat = GETINT('OutputFormat','2')
skip         = GETINT('SkipSample','1')
doFFT        = GETLOGICAL('doFFT','F')
doPSD        = GETLOGICAL('doPSD','T')
IF(doFFT.OR.doPSD) doSpec=.TRUE.
IF(doSpec) THEN
  ! two readin "modes" for spectrum averaging:
  ! 1. Prescription of number of blocks
  nBlocks=GETINT('nBlocks','1')
  ! 2. Prescription of Sampling Frequency and Blocksize
  samplingFreq=GETREAL('SamplingFreq','-999')
  cutoffFreq=GETREAL('cutoffFreq','-999.')
  IF(samplingFreq.GT.0.) THEN
    BlockSize=GETINT('BlockSize') 
  END IF
  doHanning   = GETLOGICAL('hanning','F')
  fourthDeriv = GETLOGICAL('fourthDeriv','F')
  ThirdOct    = GETLOGICAL('ThirdOct','F')
  IF (ThirdOct) THEN
    u_infPhys = GETREAL('u_inf') !velocity for re-dimensionalization of frequency
    chordPhys = GETREAL('chord')   !length for re-dimensionalization of frequency
  END IF
END IF
equiTimeSpacing=.TRUE.

!======================================================
! Readin StochInput.h5
!======================================================
CALL GET_COMMAND_ARGUMENT(2,StochFile)
CALL OpenDataFile(StochFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
!CALL ReadAttribute(File_ID,'ProjectName'  ,1,StrScalar   = ProjectName)
CALL ReadAttribute(File_ID,'nGlobalRuns'  ,1,IntScalar   = nStochSamples)
CALL ReadAttribute(File_ID,'nStochVars'   ,1,IntScalar   = nStochVars)
CALL ReadAttribute(File_ID,'polyDeg'      ,1,IntScalar   = M)

ALLOCATE(StochVarNames(nStochVars))
ALLOCATE(Distributions(nStochVars))
CALL ReadAttribute(File_ID,'StochVarNames' ,nStochVars,StrArray = StochVarNames)
CALL ReadAttribute(File_ID,'Distributions' ,nStochVars,StrArray = Distributions)

ALLOCATE(StochPoints(nStochVars,nStochSamples))
ALLOCATE(StochWeights(nStochSamples))
ALLOCATE(DistributionProps(2,nStochVars))
CALL ReadArray(ArrayName  = 'Samples',&
               Rank       = 2,&
               nVal       = (/nStochVars,nStochSamples/),&
               Offset_in  = 0,&
               Offset_dim = 1,&
               RealArray  = StochPoints)

CALL ReadArray(ArrayName  = 'Weights',&
               Rank       = 1,&
               nVal       = (/nStochSamples/),&
               Offset_in  = 0,&
               Offset_dim = 1,&
               RealArray  = StochWeights)

CALL ReadArray(ArrayName  = 'DistributionProps',&
               Rank       = 2,&
               nVal       = (/2,nStochVars/),&
               Offset_in  = 0,&
               Offset_dim = 1,&
               RealArray  = DistributionProps)

CALL CloseDataFile()

!======================================================
! Readin RP Files
!======================================================
nRPFiles = nArgs-2
ALLOCATE(RPFiles(1:nRPFiles))
DO iArg=3,nArgs
  CALL GET_COMMAND_ARGUMENT(iArg,RPFiles(iArg-2))
END DO

!======================================================
! Initialize Multiindex and polynomials
!======================================================
CALL GetHermiteCoefficients()
CALL Binom(nStochVars+M,nStochVars,nStochCoeffs)
CALL CreateMultiIndex()

END SUBROUTINE InitNisp_RP

!----------------------------------------------------------------------------------------------------------------------------------!
! Perform FFT of each RPFile and for each stochastic sample
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE PerformSampleFFT() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Globals
USE MOD_RPData                      ,ONLY: ReadRPData,AssembleRPData,FinalizeRPData
USE MOD_RPData_Vars                 
USE MOD_Nisp_RP_Vars
USE MOD_ParametersVisu
USE MOD_RPSetVisuVisu_Vars          ,ONLY: nRP_global 
USE MOD_OutputRPVisu_Vars           ,ONLY: RPData_out,nSamples_out
USE MOD_OutputRPVisu                ,ONLY: InitOutput
USE MOD_Spec                        ,ONLY: InitSpec,spec,FinalizeSpec
USE MOD_spec_Vars                        
USE MOD_RPInterpolation
USE MOD_RPInterpolation_Vars
USE MOD_EquationRP
USE MOD_RPSetVisu                   ,ONLY: FinalizeRPSet
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
! Space-separated list of input and output types. Use: (int|real|logical|...)_(in|out|inout)_dim(n)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER               :: iStochSample,iRPFile,stochOffset,iExt
CHARACTER(LEN=255)    :: InputRPFile
!===================================================================================================================================
DO iStochSample=1,nStochSamples
   ! readin RP Data from all input files
   DO iRPFile=1,nRPFiles
     InputRPFile=RPFiles(iRPFile)
     WRITE(UNIT_stdOut,'(132("="))')
     WRITE(UNIT_stdOut,'(A,I5,A,I5,A,I5)') ' PROCESSING FILE ',iRPFile,' of ',nRPFiles,' FILES'
     WRITE(UNIT_stdOut,'(A,A,A)') ' ( "',TRIM(InputRPFile),'" )'
     WRITE(UNIT_stdOut,'(132("="))')

     ! Get start index of file extension to check if it is a h5 file
     iExt=INDEX(InputRPFile,'.',BACK = .TRUE.)
     IF(InputRPFile(iExt+1:iExt+2) .NE. 'h5') &
       CALL Abort(__STAMP__,'ERROR - Invalid file extension!')

     ! Read in main attributes from given HDF5 State File
     IF(iRPFile.EQ.1) THEN
       CALL ReadRPData(InputRPFile,firstFile=.TRUE.,stochOffset=iStochSample-1)
     ELSE 
       CALL ReadRPData(InputRPFile,firstFile=.FALSE.,stochOffset=iStochSample-1)
     END IF
   END DO
   
   CALL AssembleRPData()
   IF (iStochSample .EQ. 1) THEN
     CALL InitEquationRP()
     CALL InitInterpolation()
     CALL InitSpec()
     CALL InitOutput()
     
     ALLOCATE(UTimeseries(1:nStochSamples,1:nVarVisu,nRP_global,nSamples_out))
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
   uTimeseries(iStochSample,:,:,:)= RPData_out
   
   CALL spec()
   IF (iStochSample .EQ. 1) THEN
     ALLOCATE(UFFT(1:nStochSamples,1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UMeanFFT(1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UVarFFT(1:nVarVisu,nRP_global,nSamples_spec))
     ALLOCATE(UModeFFT(1:nVarVisu,nRP_global,nSamples_spec))
     UFFT =0.
     UMeanFFT =0.
     UVarFFT =0.
     UModeFFT =0.
   END IF
     
   UFFT(iStochSample,:,:,:)=RPData_spec
   IF(iStochSample.LT.nStochSamples) THEN
     CALL FinalizeRPSet()
     CALL FinalizeRPData()
     CALL FinalizeSpec()
   END IF
END DO
time = RPTime
END SUBROUTINE PerformSampleFFT



!----------------------------------------------------------------------------------------------------------------------------------!
! Compute stochastic modes
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE ComputeModes()
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_Nisp_RP_Vars
USE MOD_EOS            ,ONLY: DefineParametersEOS,InitEOS, ConsToPrim
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataProps,ReadAttribute, ReadArray, GetDataSize
USE MOD_Basis          ,ONLY: LegendrePolynomialAndDerivative
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: jStochSample,iStochVar,l,iElem,iStochCoeff,j
REAL                                 :: evalPoly, y_out, y_out_dummy
!-----------------------------------------------------------------------------------------------------------------------------------
UMeanFFT = 0.
UVarTimeseries = 0.
UVarFFT = 0.
DO iStochCoeff=0,nStochCoeffs
  IF(iStochCoeff==0) THEN
    DO j=1, nStochSamples
      UMeanTimeseries = UMeanTimeseries + uTimeseries(j,:,:,:)*StochWeights(j)
      UMeanFFT        = UMeanFFT        + UFFT(j,:,:,:)*StochWeights(j)
    END DO
  ELSE
    y_out=0.
    UModeTimeseries = 0.
    UModeFFT = 0.
    DO jStochSample=1, nStochSamples
      evalPoly = 1.
      DO iStochVar=1, nStochVars
        IF(STRICMP(Distributions(iStochVar),"Uniform")) THEN 
          CALL LegendrePolynomialAndDerivative(MultiIndex(iStochCoeff,iStochVar), &
              (2*StochPoints(iStochVar,jStochSample)-(DistributionProps(2,iStochVar)+DistributionProps(1,iStochVar)))&
              /(DistributionProps(2,iStochVar)-DistributionProps(1,iStochVar)),y_out,y_out_dummy)
        ELSE !Hermite
          CALL EvaluateHermitePoly(MultiIndex(iStochCoeff,iStochVar),&
              (StochPoints(iStochVar,jStochSample)-DistributionProps(1,iStochVar))/DistributionProps(2,iStochVar),y_out)
        END IF
        evalPoly = evalPoly*y_out
      END DO
      UModeTimeseries = UModeTimeseries+ uTimeseries(jStochSample,:,:,:)*evalPoly*StochWeights(jStochSample)
      UModeFFT = UModeFFT+ UFFT(jStochSample,:,:,:)*evalPoly*StochWeights(jStochSample)
    END DO  
  END IF
  UVarTimeseries = UVarTimeseries + UModeTimeseries*UModeTimeseries
  UVarFFT = UVarFFT + UModeFFT*UModeFFT
END DO
END SUBROUTINE ComputeModes

!----------------------------------------------------------------------------------------------------------------------------------!
! Finalize all parameters of NISP RP tool
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeNisp_RP() 
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Nisp_RP_Vars 
USE MOD_Spec            ,ONLY: FinalizeSpec
USE MOD_EquationRP      ,ONLY: FinalizeEquationRP
USE MOD_RPInterpolation ,ONLY: FinalizeInterpolation
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL FinalizeInterpolation()
CALL FinalizeEquationRP()
CALL FinalizeSpec()

SDEALLOCATE(UMeanTimeseries)
SDEALLOCATE(UVarTimeseries)
SDEALLOCATE(UModeTimeseries) 
SDEALLOCATE(UFFT)
SDEALLOCATE(UMeanFFT)
SDEALLOCATE(UVarFFT)
SDEALLOCATE(UModeFFT)
SDEALLOCATE(time)
SDEALLOCATE(UTimeseries)
SDEALLOCATE(StochVarNames)
SDEALLOCATE(Distributions)
SDEALLOCATE(StochPoints)
SDEALLOCATE(StochWeights)
SDEALLOCATE(DistributionProps)
SDEALLOCATE(RPFiles)
SDEALLOCATE(HermiteCoeffs)
SDEALLOCATE(MultiIndex)
END SUBROUTINE FinalizeNisp_RP


!===================================================================================================================================
!> Create Multi-Index to evaluate Polynomials
!> M = highest polynomial degree |i|=M, i=(i1,...,i_stochDim)
!> nStochCoeffs = #number of stochdim-dimensional chaos polynomials of order <=M
!> Algorithm in LeMaitre, p.524
!===================================================================================================================================
SUBROUTINE CreateMultiIndex()
!MODULES
USE MOD_Nisp_RP_Vars, ONLY: MultiIndex, nStochVars, M, nStochCoeffs
IMPLICIT NONE!
!----------------------------------------------------------------------------------------------------------------------------------!
! LOCAL VARIABLES
INTEGER           :: i,j,k,p,q,r 
INTEGER           :: p2(1:M,1:nStochVars)
!----------------------------------------------------------------------------------------------------------------------------------!
ALLOCATE(MultiIndex(0:nStochCoeffs,1:nStochVars))
MultiIndex = 0.
DO i=1,nStochVars
  MultiIndex(i,i) = 1
END DO
p2= 0.
q= nStochVars
p2(1,:)=1
DO k=2,M
  r=q
  DO i=1,nStochVars
    DO p=i,nStochVars
      p2(k,i)= p2(k,i) + p2(k-1,p)
    END DO
  END DO
  DO j=1,nStochVars
    DO p= r- p2(k,j)+1, r
      q=q+1
      MultiIndex(q,:)=MultiIndex(p,:)
      MultiIndex(q,j)=MultiIndex(q,j)+1
    END DO
  END DO
END DO
END SUBROUTINE CreateMultiIndex

!===================================================================================================================================
!> Get the hermite coeffficients of the Hermite polynomials
!===================================================================================================================================
SUBROUTINE GetHermiteCoefficients()
! MODULES
USE MOD_Nisp_RP_Vars, ONLY:HermiteCoeffs,M
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p
!==================================================================================================================================
ALLOCATE(HermiteCoeffs(0:M,0:M))
HermiteCoeffs(:,:)=0.
HermiteCoeffs(0,0)=1.
HermiteCoeffs(1,1)=1.
DO p=2,M
  HermiteCoeffs(1:p,p)   = HermiteCoeffs(0:p-1,p-1)/SQRT(REAL(p))
  HermiteCoeffs(0:p-2,p) = HermiteCoeffs(0:p-2,p)-HermiteCoeffs(0:p-2,p-2)*SQRT(REAL(p-1)/REAL(p))
END DO
END SUBROUTINE GetHermiteCoefficients


!===================================================================================================================================
!> evaluate Hermite polynomial of degree pStoch at xEval
!===================================================================================================================================
SUBROUTINE EvaluateHermitePoly(pStoch,xEval,y_out)
! MODULES
USE MOD_Nisp_RP_Vars, ONLY:HermiteCoeffs 
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
  y_out=y_out+HermiteCoeffs(p,pStoch)*xEval**p  
END DO
END SUBROUTINE EvaluateHermitePoly


!===================================================================================================================================
!> Compute the binomial coefficient of n over k where n=M+nStochVars and k=M
!===================================================================================================================================
SUBROUTINE Binom(n, k, resu)
! MODULES
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES 
INTEGER , INTENT(IN)   :: n,k
INTEGER , INTENT(OUT)  :: resu
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: i,j
!-----------------------------------------------------------------------------------------------------------------------------------
j=k
resu=1
IF(k.GT.(n-k)) THEN
  j = n - k
END IF
DO i=0,j-1
  resu = resu*(n - i)
  resu = resu/(i + 1)
END DO
END SUBROUTINE Binom

END MODULE MOD_NISP_RP
