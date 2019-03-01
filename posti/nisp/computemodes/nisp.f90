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
!===================================================================================================================================
! Add comments please!
!===================================================================================================================================


PROGRAM NISP
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Nisp_Vars
USE MOD_Nisp_Output
USE MOD_ReadInTools
USE MOD_StringTools             ,ONLY: STRICMP,GetFileExtension
USE MOD_MPI                     ,ONLY: DefineParametersMPI,InitMPI
USE MOD_Interpolation           ,ONLY: DefineParametersInterpolation,InitInterpolation,FinalizeInterpolation
USE MOD_Interpolation_Vars      ,ONLY: xGP,wBary
USE MOD_IO_HDF5                 ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI                     ,ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_IO_HDF5                 ,ONLY: AddToFieldData
USE MOD_Mesh                    ,ONLY: DefineParametersMesh,InitMesh
USE MOD_Mesh_Vars               ,ONLY: nGlobalElems
USE MOD_Exactfunc               ,ONLY: ExactFunc
USE MOD_Analyze                 ,ONLY: InitAnalyzeBasis
USE MOD_Output_Vars             ,ONLY: ProjectName
USE MOD_EOS                     ,ONLY: ConsToPrim
USE MOD_EOS_Vars                ,ONLY: KappaM1,R, Kappa
USE MOD_IO_HDF5                 ,ONLY: File_ID,nDims,HSize
USE MOD_HDF5_Input              ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataProps,ReadAttribute, ReadArray, GetDataSize
USE MOD_HDF5_Output             ,ONLY: WriteAttribute,WriteArray,GenerateFileSkeleton,WriteHeader
USE MOD_ReadInTools             ,ONLY: ExtractParameterFile
USE MOD_Output_Vars             ,ONLY: UserBlockTmpFile,userblock_total_len
USE MOD_Output                  ,ONLY: insert_userblock
USE ISO_C_BINDING               ,ONLY: C_NULL_CHAR
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=2)                     :: integer_string
CHARACTER(LEN=7)                     :: disfunc= "DisFunc"
CHARACTER(LEN=7)                     :: disprop= "DisProp"
CHARACTER(LEN=9)                     :: distribution_string
INTEGER                              :: i, res, nVar_Field
INTEGER                              :: line_no
LOGICAL                              :: FieldDataFound, create
LOGICAL                              :: userblockFound
CHARACTER(LEN=255)                   :: prmfile=".parameter.ini"
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()



SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(1("**********************************"))')
SWRITE(UNIT_stdOut,'(1("**     START COMPUTING MODES    **"))')
SWRITE(UNIT_stdOut,'(1("**********************************"))')
SWRITE(UNIT_stdOut,'(132("="))')

CALL prms%SetSection("NISP Setup")
CALL prms%SetSection("Stochastic Input")
CALL prms%CreateIntOption( "NUncertainties"             , "Number of stochastic input variables", multiple=.TRUE.)
CALL prms%CreateIntOption( "NStoch"             , "max. polynomial chaos degree",multiple=.TRUE.)
CALL prms%CreateIntOption(      "TotalSamples"             , "Number of samples which need to be computedl",multiple=.TRUE.)
CALL prms%SetSection("Equation of State")
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent",multiple=.TRUE.)
CALL prms%CreateRealOption(     'R',            "Specific gas constant",multiple=.TRUE.)
CALL prms%read_options(Args(1))
SWRITE(UNIT_stdOut,'(132("="))')

IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: parameter_nisp.ini statefile_last_sample')
END IF
CALL InitIOHDF5()
! CALL InitEos()
nStochVar  = GETINT('NUncertainties')
N0         = GETINT('NStoch')
totSamples = GETINT('TotalSamples','-1')
IF(totSamples.LT.0) THEN
  totSamples  = N0+1
END IF
Kappa    =GETREAL('kappa','1.4')
KappaM1  =Kappa-1.
R=GETREAL('R','287.058')

ALLOCATE(distributions(1:nStochVar))
ALLOCATE(weights(1:totSamples))
ALLOCATE(quadpoints(1:totSamples, 1:nStochVar))
ALLOCATE(quadpointsTransformed(1:totSamples, 1:nStochVar))
ALLOCATE(quadpointsInv(1:totSamples, 1:nStochVar))
! Open .h5 on sample n to get MeshFile and necessary attributes
StateFileName=TRIM(Args(2))
WRITE(sample_string,'(I0)') totSamples
StateFile= 'sample_'//TRIM(sample_string)//'/'//TRIM(StateFileName)

DO i=1,nStochVar
  write (integer_string,'(I0)') i
  distribution_string = disfunc//integer_string
  CALL prms%CreateIntOption(      distribution_string             , "1: uniform distribution  2: Normal distributuion")
END DO

CALL prms%read_options(Args(1))

DO i=1,nStochVar
  write (integer_string,'(I0)') i
  distribution_string = disfunc//integer_string
  distributions(i)=  GETINT(distribution_string)
END DO


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

DO i=1,nStochVar
  IF(distributions(i)==2) THEN
    CALL GetHermiteCoefficients()
    EXIT
  END IF
END DO

CALL binom(nStochVar+N0,nStochVar, nPoly)

CALL CreateMultiIndex()


!======================================================
! Open .h5 on sample n to get MeshFile and necessary attributes
!======================================================
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ExtractParameterFile(StateFile,TRIM(prmfile),userblockFound)
 !prepare userblock file
CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(prmfile)//C_NULL_CHAR)
INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)

CALL GetDataProps(nVar,NNew,nElemsNew,NodeType)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=Mesh)
CALL ReadAttribute(File_ID,'Project_Name',1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_State)

IF(.NOT. ALLOCATED(VarNames)) ALLOCATE(VarNames(nVar))
CALL ReadAttribute(File_ID,'VarNames',nVar,StrArray=VarNames)
CALL DatasetExists(File_ID, 'FieldData', FieldDataFound)
IF (FieldDataFound) THEN
  ! get size of FieldData array
  CALL GetDataSize(File_ID,'FieldData',nDims,HSize)
  nVar_Field=INT(HSize(1),4)
  IF(.NOT. ALLOCATED(VarNamesField)) ALLOCATE(VarNamesField(nVar_Field))
  CALL ReadAttribute(File_ID,'VarNamesAddField',nVar_Field,StrArray=VarNamesField)
END IF

CALL CloseDataFile()



PP_N= NNew
nGlobalElems = nElemsNew
ALLOCATE(UMean(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew), UVar(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(U(nVar,0:NNew,0:NNew,0:NNew,nElemsNew),UMode(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(UPrim(nVar+1,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
UMean= 0.
UVar = 0.
UMode = 0.
UPrim =0.
!----------------------------------------------------------------------------------------------------------------------------------!
!EXPECTATION
!----------------------------------------------------------------------------------------------------------------------------------!
CALL ComputeMode(0)
!----------------------------------------------------------------------------------------------------------------------------------!
!VARIANCE
!----------------------------------------------------------------------------------------------------------------------------------!
DO i=1, nPoly-1
  CALL ComputeMode(i)
  UVar = UVar + UMode*UMode
END DO

SWRITE(UNIT_stdOut,'(A)') ' WRITING  MEAN and VARIANCE...'
CALL WriteMeanAndVarianceToHDF5()

SDEALLOCATE(UMean)
SDEALLOCATE(UVar)
SDEALLOCATE(U)
SDEALLOCATE(UMode)
SDEALLOCATE(UPrim)
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' NISP TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM NISP

SUBROUTINE ComputeMode(i)
USE MOD_Nisp_Vars
USE MOD_EOS                     ,ONLY: DefineParametersEOS,InitEOS, ConsToPrim
USE MOD_HDF5_Input     ,ONLY: OpenDataFile,CloseDataFile,DatasetExists,GetDataProps,ReadAttribute, ReadArray, GetDataSize
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER, INTENT(IN) :: i
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER           :: j,k,l,m,iElem
REAL              :: evalPoly, y_out, y_out_dummy
!-----------------------------------------------------------------------------------------------------------------------------------
IF(i==0) THEN
  y_out=0.
  DO j=1, totSamples
    WRITE(sample_string,'(I0)') j
    StateFile= 'sample_'//TRIM(sample_string)//'/'//TRIM(StateFileName)
    ! Open the data file
    CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    ! Read the DG solution and store in UNew
    CALL ReadArray('DG_Solution',5,&
               (/nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),0,5,RealArray=U)
    DO iElem=1,nElemsNew
      DO k=0,NNew; DO l=0,NNew; DO m=0,NNew
        CALL ConsToPrim(UPrim(:,m,l,k,iElem),U(:,m,l,k,iElem))
      END DO; END DO; END DO
    END DO
    UMean(1:nVar,:,:,:,:) = UMean(1:nVar,:,:,:,:) + U(:,:,:,:,:)*weights(j)
    UMean(nVar+1:,:,:,:,:) = UMean(nVar+1:,:,:,:,:) + UPrim(2:,:,:,:,:)*weights(j)
  END DO
ELSE
  UMode = 0.
  DO j=1, totSamples
    evalPoly = 1.
    DO k=1, nStochVar
      IF(distributions(k)==1) THEN !Legendre on [-0.5,0.5]
        CALL LegendrePolynomialAndDerivative(MultiIndex(i,k), (quadpoints(j,k)-0.5),.TRUE.,y_out,y_out_dummy)
      ELSE !Hermite
        CALL EvaluateHermitePoly(MultiIndex(i,k),quadpoints(j,k),y_out)
      END IF
      evalPoly = evalPoly*y_out
    END DO
    WRITE(sample_string,'(I0)') j
    StateFile= 'sample_'//TRIM(sample_string)//'/'//TRIM(StateFileName)
    ! Open the data file
    CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
    ! Read the DG solution and store in UNew
    CALL ReadArray('DG_Solution',5,&
               (/nVar,NNew+1,NNew+1,NNew+1,nElemsNew/),0,5,RealArray=U)
    DO iElem=1,nElemsNew
      DO k=0,NNew; DO l=0,NNew; DO m=0,NNew
        CALL ConsToPrim(UPrim(:,m,l,k,iElem),U(:,m,l,k,iElem))
      END DO; END DO; END DO
    END DO
    UMode(1:nVar,:,:,:,:) = UMode(1:nVar,:,:,:,:) + U*evalPoly*weights(j)
    UMode(nVar+1:,:,:,:,:)  = UMode(nVar+1:,:,:,:,:)  + UPrim(2:,:,:,:,:)*evalPoly*weights(j)
  END DO
END IF


END SUBROUTINE

SUBROUTINE CreateMultiIndex()
! INPUT / OUTPUT VARIABLES
!N0 = highest polynomial degree |i|=N0, i=(i1,...,i_stochDim)
!nPoly = #number of stochdim-dimensional chaos polynomials of order <=N0
!Algorithm in LeMaitre, p.524
!MODULES
USE MOD_Nisp_Vars, ONLY: MultiIndex, nStochVar, N0, nPoly
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
USE MOD_Nisp_Vars, ONLY:HermiteCoeff,N0
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
USE MOD_Nisp_Vars, ONLY:HermiteCoeff
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
