!=================================================================================================================================
! Copyright (c) 2010-2022  Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024  Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
#include "eos.h"

!==================================================================================================================================
!> Contains analyze routines specific to the linear scalar advection equation
!==================================================================================================================================
MODULE MOD_AnalyzeEquation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitAnalyzeEquation
  MODULE PROCEDURE InitAnalyzeEquation
END INTERFACE

INTERFACE AnalyzeEquation
  MODULE PROCEDURE AnalyzeEquation
END INTERFACE

INTERFACE FinalizeAnalyzeEquation
  MODULE PROCEDURE FinalizeAnalyzeEquation
END INTERFACE


PUBLIC:: AnalyzeEquation, InitAnalyzeEquation, FinalizeAnalyzeEquation
PUBLIC:: DefineParametersAnalyzeEquation
!==================================================================================================================================

CONTAINS
   
!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("AnalyzeEquation")
CALL prms%CreateLogicalOption('CalcAcousticEnergy'   , "Set true to compute global acoustic energy norm"         , '.FALSE.')
CALL prms%CreateLogicalOption('WriteAcousticEnergy'  , "Set true to write acoustic energy to file"             , '.TRUE.')
END SUBROUTINE DefineParametersAnalyzeEquation


!==================================================================================================================================
!> Initializes variables necessary for analyse subroutines
!==================================================================================================================================
SUBROUTINE InitAnalyzeEquation()
! MODULES
USE MOD_Globals
USE MOD_AnalyzeEquation_Vars
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_Output,             ONLY: InitOutputToFile
USE MOD_Output_Vars,        ONLY: ProjectName
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
doCalcAcousticEnergy    =GETLOGICAL('CalcAcousticEnergy'   ,'.FALSE.')
doWriteAcousticEnergy   =GETLOGICAL('WriteAcousticEnergy'  ,'.TRUE.')
! Initialize eval routines
IF(MPIRoot)THEN
  IF(doCalcAcousticEnergy.AND.doWriteAcousticEnergy)THEN
    FileName_AcousticEnergy  = TRIM(ProjectName)//'_acE'
    CALL InitOutputToFile(FileName_AcousticEnergy,'AcousticEnergy',1,(/'AcousticEnergy'/))
  END IF
END IF
END SUBROUTINE InitAnalyzeEquation


!==================================================================================================================================
!> Equation specific analyze routine
!==================================================================================================================================
SUBROUTINE AnalyzeEquation(Time)
! MODULES
USE MOD_Globals
USE MOD_AnalyzeEquation_Vars
USE MOD_Output,             ONLY: OutputToFile
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: Time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=40)               :: formatStr
REAL :: meanE
!==================================================================================================================================
IF(doCalcAcousticEnergy)  CALL CalcAcousticEnergy(meanE)

IF(MPIRoot.AND.doCalcAcousticEnergy)THEN
  IF (doWriteAcousticEnergy) &
    CALL OutputToFile(FileName_AcousticEnergy,(/Time/),(/1,1/),(/meanE/))
  WRITE(formatStr,'(A,I2,A)')'(A14,',1,'ES18.9)'
  WRITE(UNIT_StdOut,formatStr)' Mean acoustic energy : ',meanE
END IF

END SUBROUTINE AnalyzeEquation



!==================================================================================================================================
!> Calculates mean acoustic energy in the whole domain. Using the definition E=1/2\rho_0 \vec{u}'^2 + 1/2\frac{p'^2}{\rho_0 c_0^2}
!> For resting base flows without sources, this quantity is conserved.
!==================================================================================================================================
SUBROUTINE CalcAcousticEnergy(meanE)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Analyze_Vars,       ONLY: wGPVol,Vol
USE MOD_Mesh_Vars,          ONLY: sJ,nElems
USE MOD_DG_Vars,            ONLY: U
USE MOD_Baseflow_Vars,      ONLY: BaseState
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)                :: meanE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: IntegrationWeight
INTEGER                         :: iElem,i,j,k
!==================================================================================================================================
meanE=0.
DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      IntegrationWeight=wGPVol(i,j,k)/sJ(i,j,k,iElem,0)
      meanE                =meanE +IntegrationWeight*(0.5*BaseState(DENS)*SUM(U(VELV,i,j,k,iElem)**2) &
                                             +0.5/(BaseState(DENS)*BaseState(SOSP)**2)*U(PRES,i,j,k,iElem)**2)
    END DO; END DO; END DO !i,j,k
END DO ! iElem

#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,meanE,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
ELSE
  CALL MPI_REDUCE(meanE           ,0  ,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif

meanE=meanE/Vol
END SUBROUTINE CalcAcousticEnergy


!==================================================================================================================================
!> Finalizes variables necessary for analyze subroutines
!==================================================================================================================================
SUBROUTINE FinalizeAnalyzeEquation()
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeAnalyzeEquation

END MODULE MOD_AnalyzeEquation
