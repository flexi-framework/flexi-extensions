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

MODULE MOD_MC

INTERFACE DefineParametersMC
  MODULE PROCEDURE DefineParametersMC
END INTERFACE

INTERFACE InitMC
  MODULE PROCEDURE InitMC
END INTERFACE

INTERFACE ComputeMCEstimator
  MODULE PROCEDURE ComputeMCEstimator
END INTERFACE

INTERFACE FinalizeMC
  MODULE PROCEDURE FinalizeMC
END INTERFACE

PUBLIC::DefineParametersMC,InitMC,ComputeMCEstimator,FinalizeMC
CONTAINS
!===================================================================================================================================
!> Define parameters of Posti MC RP tool
!===================================================================================================================================
SUBROUTINE DefineParametersMC()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_ReadInTools
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL prms%SetSection('MC Parameters')
CALL prms%CreateRealArrayOption("RefState"           , " RefState, for derived quantities")
END SUBROUTINE DefineParametersMC

!===================================================================================================================================
!> Read parameters of Posti MC RP tool
!===================================================================================================================================
SUBROUTINE InitMC()
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Readintools
USE MOD_MC_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input        ,ONLY: OpenDataFile,CloseDataFile,GetDataProps,ReadAttribute, ReadArray, ReadAttributeBatchScalar
USE MOD_HDF5_Input        ,ONLY: ReadAttribute,ReadArray,ReadAttributeBatchScalar,GetDataSize
USE MOD_Output            ,ONLY: insert_userblock
USE MOD_Output_Vars       ,ONLY: UserBlockTmpFile,userblock_total_len
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
LOGICAL                 :: userblockFound
CHARACTER(LEN=255)      :: prmfile=".parameter.ini"
!===================================================================================================================================
!======================================================
! Open .h5 on sample n to get MeshFile and necessary attributes
!======================================================
CALL GET_COMMAND_ARGUMENT(2,StateFileName)
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
CALL ReadAttribute(File_ID,'Project_Name',    1,StrScalar=ProjectName)
CALL ReadAttribute(File_ID,'MeshFile',    1,StrScalar=Mesh)
CALL ReadAttribute(File_ID,'Time',1,RealScalar=Time_State)

CALL CloseDataFile()

addVars=1
ALLOCATE(UMean(2*nVar+addVars,0:NNew,0:NNew,0:NNew,nElemsNew), UVar(2*nVar+addVars,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(UMeanSquare(2*nVar+addVars,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(U(nVar,0:NNew,0:NNew,0:NNew,nElemsNew))
ALLOCATE(UPrim(nVar+1,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
ALLOCATE(UAddVars(addVars,0:NNew,  0:NNew,  0:NNew,  nElemsNew))
UMean = 0.
UVar= 0.
UMeanSquare= 0.
U = 0.
UPrim = 0.
UAddVars = 0.
nSamples=0
RefState = GETREALARRAY('RefState',nVar)
END SUBROUTINE InitMC

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadStateFile(StateFileLocal,DataSetName,offset)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray
USE MOD_MC_Vars
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

!===================================================================================================================================
!> Open a state file, read the old state and store the information later needed to write a new state.
!===================================================================================================================================
SUBROUTINE ReadNumberOfSamples(StateFileLocal,DataSetName)
! MODULES                                                                                                                          !
USE MOD_HDF5_Input,    ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataSize
USE MOD_IO_HDF5
USE MOD_MC_Vars
USE MOD_Globals
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
CHARACTER(LEN=255),INTENT(IN) :: StateFileLocal !< State file to be read
CHARACTER(LEN=255),INTENT(IN) :: DataSetName
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL OpenDataFile(StateFileLocal,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
SWRITE(UNIT_stdOut,'(a,a,a)',ADVANCE='NO')' Read Data from HDF5 FILE "',TRIM(StateFileLocal),'" ... \n'
CALL GetDataSize(File_ID,DataSetName,nDims,HSize)
iFileSamples=INT(HSize(nDims))
DEALLOCATE(HSize)
CALL CloseDataFile()
END SUBROUTINE ReadNumberOfSamples

!----------------------------------------------------------------------------------------------------------------------------------!
! Compute stochastic modes
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE ComputeMCEstimator()
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_MC_Vars
USE MOD_StringTools    ,ONLY: STRICMP
USE MOD_EOS            ,ONLY: DefineParametersEOS,InitEOS, ConsToPrim
USE MOD_Basis          ,ONLY: LegendrePolynomialAndDerivative
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: jStochSample
CHARACTER(LEN=255)      :: DataSetName="DG_Solution"
INTEGER                 :: i,k,l,m,iElem
!-----------------------------------------------------------------------------------------------------------------------------------
DO i=1,nStateFiles
  CALL GET_COMMAND_ARGUMENT(i+1,StateFileName)
  CALL ReadNumberOfSamples(StateFileName,DataSetName)
  nSamples=nSamples+iFileSamples
  DO jStochSample=1,iFileSamples
    CALL ReadStateFile(StateFileName,DataSetName,jStochSample-1)
    DO iElem=1,nElemsNew
      DO k=0,NNew; DO l=0,NNew; DO m=0,NNEW
        CALL ConsToPrim(UPrim(:,m,l,k,iElem),U(:,m,l,k,iElem))
      END DO; END DO; END DO
    END DO
    UAddVars(1,:,:,:,:) = (UPrim(5,:,:,:,:)-RefState(5))/(0.5*RefState(1)*NORM2(RefState(2:4))*NORM2(RefState(2:4)))
    UMean(1:nVar,:,:,:,:)                 = UMean(1:nVar,:,:,:,:)                  + U
    UMean(nVar+1:2*nVar,:,:,:,:)          = UMean(nVar+1:2*nVar,:,:,:,:)           + UPrim(2:nVar+1,:,:,:,:) 
    UMean(2*nVar+1:2*nVar+addVars,:,:,:,:)= UMean(2*nVar+1:2*nVar+addVars,:,:,:,:) + UAddVars 
    UMeanSquare(1:nVar,:,:,:,:)                 = UMeanSquare(1:nVar,:,:,:,:)                  + U**2.
    UMeanSquare(nVar+1:2*nVar,:,:,:,:)          = UMeanSquare(nVar+1:2*nVar,:,:,:,:)           + UPrim(2:nVar+1,:,:,:,:)**2. 
    UMeanSquare(2*nVar+1:2*nVar+addVars,:,:,:,:)= UMeanSquare(2*nVar+1:2*nVar+addVars,:,:,:,:) + UAddVars**2.0
  END DO
END DO 

UMean = UMean/nSamples
UMeanSquare = UMeanSquare/nSamples
UVar = UMeanSquare - UMean**2.0
END SUBROUTINE ComputeMCEstimator

!----------------------------------------------------------------------------------------------------------------------------------!
! Finalize all parameters of NISP tool
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE FinalizeMC()
!----------------------------------------------------------------------------------------------------------------------------------!
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_MC_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SDEALLOCATE(UMean)
SDEALLOCATE(UVar)
SDEALLOCATE(U)
SDEALLOCATE(UPrim)
SDEALLOCATE(UAddVars)
SDEALLOCATE(UMeanSquare)
SDEALLOCATE(StochVarNames)
END SUBROUTINE FinalizeMC

END MODULE MOD_MC
