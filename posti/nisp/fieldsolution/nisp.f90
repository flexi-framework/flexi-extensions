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

MODULE MOD_NISP

INTERFACE DefineParametersNisp
  MODULE PROCEDURE DefineParametersNisp
END INTERFACE

INTERFACE InitNisp
  MODULE PROCEDURE InitNisp
END INTERFACE

INTERFACE ComputeModes
  MODULE PROCEDURE ComputeModes
END INTERFACE

INTERFACE FinalizeNisp
  MODULE PROCEDURE FinalizeNisp
END INTERFACE

PUBLIC::DefineParametersNisp,InitNisp,ComputeModes,FinalizeNisp
CONTAINS
!===================================================================================================================================
!> Define parameters of Posti Nisp RP tool
!===================================================================================================================================
SUBROUTINE DefineParametersNisp()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ReadInTools
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%SetSection('Nisp Parameters')
CALL prms%CreateRealOption   (  'kappa'              , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'Pr'                 , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'R'                  , "TODO",multiple=.TRUE.)
CALL prms%CreateRealOption   (  'mu0'                , "TODO",multiple=.TRUE.)
END SUBROUTINE DefineParametersNisp

!===================================================================================================================================
!> Read parameters of Posti Nisp RP tool
!===================================================================================================================================
SUBROUTINE InitNisp()
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Readintools
USE MOD_Nisp_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input        ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute, ReadArray, ReadAttributeBatchScalar
USE MOD_HDF5_Input        ,ONLY: ReadAttribute,ReadArray,ReadAttributeBatchScalar,GetDataSize
USE MOD_Output            ,ONLY: insert_userblock
USE MOD_Output_Vars       ,ONLY: UserBlockTmpFile,userblock_total_len
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255)      :: StochFile
INTEGER                 :: iArg,nLevelVarsStr
LOGICAL                 :: userblockFound
CHARACTER(LEN=255)      :: prmfile=".parameter.ini"
!===================================================================================================================================
!======================================================
! Readin StochInput.h5
!======================================================
CALL GET_COMMAND_ARGUMENT(2,StochFile)
CALL OpenDataFile(StochFile,create=.FALSE.,single=.TRUE.,readOnly=.TRUE.)
!CALL ReadAttribute(File_ID,'ProjectName'  ,1,StrScalar   = ProjectName)
CALL ReadAttribute(File_ID,'nGlobalRuns'  ,1,IntScalar   = nStochSamples)
CALL ReadAttribute(File_ID,'nStochVars'   ,1,IntScalar   = nStochVars)
CALL ReadAttribute(File_ID,'polyDeg'      ,1,IntScalar   = M)
CALL ReadAttributeBatchScalar('ProjectName',StrScalar = ProjectName)

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
! Initialize Multiindex and polynomials
!======================================================
CALL GetHermiteCoefficients()
CALL Binom(nStochVars+M,nStochVars,nStochCoeffs)
nStochCoeffs = nStochCoeffs-1
CALL CreateMultiIndex()
!======================================================
! Open .h5 on sample n to get MeshFile and necessary attributes
!======================================================
CALL GET_COMMAND_ARGUMENT(3,StateFileName)
CALL OpenDataFile(StateFileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ExtractParameterFile(StateFileName,TRIM(prmfile),userblockFound)
#if FV_ENABLED
CALL GetDataSize(File_ID,'ElemData',nDims,HSize)
CHECKSAFEINT(MAXVAL(HSize),4)
nVal_ElemData=  INT(HSize)
DEALLOCATE(HSize)
ALLOCATE(ElemData(nVal_ElemData(1),nVal_ElemData(2),1))
ALLOCATE(VarNamesAdditional(nVal_ElemData(1)))
CALL ReadArray('ElemData',3,(/nVal_ElemData(1),nVal_ElemData(2),1/),0,1,RealArray=ElemData)
CALL ReadAttribute(File_ID,'VarNamesAdd',nVal_ElemData(1),StrArray=VarNamesAdditional)
#endif
 !prepare userblock file
CALL insert_userblock(TRIM(UserBlockTmpFile)//C_NULL_CHAR,TRIM(prmfile)//C_NULL_CHAR)
INQUIRE(FILE=TRIM(UserBlockTmpFile),SIZE=userblock_total_len)

CALL GetDataProps(nVar,NNew,nElemsNew,NodeType)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=Mesh)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_State)

CALL CloseDataFile()
ALLOCATE(UMean(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew), UVar(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(UMode(2*nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(U(nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(UPrim(nVar+1,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(USamplePrim(1:nStochSamples,1:nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(USampleCons(1:nStochSamples,1:nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
UMean = 0.
UVar= 0.
UMode= 0.
U = 0.
UPrim = 0.
USamplePrim = 0.
USampleCons = 0.
END SUBROUTINE InitNisp

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadStateFile(StateFileLocal,DataSetName,offset)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_IO_HDF5,       ONLY: File_ID
USE MOD_Nisp_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: StateFileLocal !< State file to be read
CHARACTER(LEN=255),INTENT(IN) :: DataSetName
INTEGER,INTENT(IN)            :: offset
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(StateFileLocal,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadArray(TRIM(DataSetName),6,&
               (/nVar,NNew+1,NNew+1,NNew+1,nElemsNew,1/),offset,6,RealArray=U)
CALL CloseDataFile()
END SUBROUTINE ReadStateFile

!----------------------------------------------------------------------------------------------------------------------------------!
! Compute stochastic modes
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE ComputeModes()
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Nisp_Vars
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_EOS            ,ONLY: DefineParametersEOS,InitEOS, ConsToPrim
USE MOD_Basis          ,ONLY: LegendrePolynomialAndDerivative
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: iStochCoeff,jStochSample,iStochVar
REAL                    :: evalPoly, y_out, y_out_dummy
CHARACTER(LEN=255)      :: DataSetName="DG_Solution"
INTEGER                 :: iStochSample,k,l,p,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
DO jStochSample=1,nStochSamples
  CALL ReadStateFile(StateFileName,DataSetName,jStochSample-1)
  USampleCons(jStochSample,:,:,:,:,:) = U
  DO iElem=1,nElemsNew
    DO k=0,NNew; DO l=0,NNew; DO p=0,NNew
      CALL ConsToPrim(UPrim(:,m,l,k,iElem),U(:,m,l,k,iElem))
    END DO; END DO; END DO
  END DO
  USamplePrim(jStochSample,:,:,:,:,:) = UPrim(2:nVar,:,:,:,:)
END DO
DO iStochCoeff=0,nStochCoeffs
  IF(iStochCoeff==0) THEN
    DO jStochSample=1, nStochSamples
      UMean(1:nVar,:,:,:,:)        = UMean(1:nVar,:,:,:,:)  + USampleCons(jStochSample,:,:,:,:,:)*StochWeights(jStochSample)
      UMean(nVar+1:2*nVar,:,:,:,:) = UMean(nVar+1:2*nVar,:,:,:,:) + USamplePrim(jStochSample,2:nVar,:,:,:,:)*StochWeights(jStochSample)
    END DO
  ELSE
    y_out=0.
    UMode =0.
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
      UMode(1:nVar,:,:,:,:)       = UMode(1:nVar,:,:,:,:) + USampleCons(jStochSample,:,:,:,:,:)*evalPoly*StochWeights(jStochSample)
      UMode(nVar+1:nVar,:,:,:,:)  = UMode(nVar+1:nVar,:,:,:,:) + USamplePrim(jStochSample,2:nVar,:,:,:,:)*evalPoly*StochWeights(jStochSample)
    END DO
   UVar = UVar + UMode*UMode
  END IF
END DO
END SUBROUTINE ComputeModes

!----------------------------------------------------------------------------------------------------------------------------------!
! Finalize all parameters of NISP tool
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeNisp()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_Nisp_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(UMean)
SDEALLOCATE(UVar)
SDEALLOCATE(U)
SDEALLOCATE(UMode)
SDEALLOCATE(USamplePrim)
SDEALLOCATE(USampleCons)
SDEALLOCATE(StochVarNames)
SDEALLOCATE(Distributions)
SDEALLOCATE(StochPoints)
SDEALLOCATE(StochWeights)
SDEALLOCATE(DistributionProps)
SDEALLOCATE(HermiteCoeffs)
SDEALLOCATE(MultiIndex)
END SUBROUTINE FinalizeNisp

!===================================================================================================================================
!> Create Multi-Index to evaluate Polynomials
!> M = highest polynomial degree |i|=M, i=(i1,...,i_stochDim)
!> nStochCoeffs = #number of stochdim-dimensional chaos polynomials of order <=M
!> Algorithm in LeMaitre, p.524
!===================================================================================================================================
SUBROUTINE CreateMultiIndex()
!MODULES
USE MOD_Nisp_Vars, ONLY: MultiIndex, nStochVars, M, nStochCoeffs
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
USE MOD_Nisp_Vars, ONLY:HermiteCoeffs,M
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
USE MOD_Nisp_Vars, ONLY:HermiteCoeffs
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

END MODULE MOD_NISP