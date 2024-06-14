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
!> This module contains the routines to calculate the equation system specific allowable timestep.
!==================================================================================================================================
MODULE MOD_CalcTimeStep
! MODULES
USE MOD_TimeDisc_Vars ,ONLY:TimeStep_global
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitCalctimestep
  MODULE PROCEDURE InitCalctimestep
END INTERFACE

INTERFACE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP
  MODULE PROCEDURE CALCTIMESTEP_dummy
END INTERFACE

INTERFACE FinalizeCalctimestep
  MODULE PROCEDURE FinalizeCalctimestep
END INTERFACE


PUBLIC :: InitCalctimestep,CALCTIMESTEP,FinalizeCalctimestep
!==================================================================================================================================
REAL,ALLOCATABLE :: MetricsAdv(:,:,:,:,:,:)  !< support variable: NORM2(Metricsfgh)/J

CONTAINS

!==================================================================================================================================
!> Precompute some metric support variables
!==================================================================================================================================
SUBROUTINE InitCalctimestep()
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
#if PARABOLIC
USE MOD_EOS_Vars ,ONLY:KappasPr
#endif /*PARABOLIC*/
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem,FVE,errType
!==================================================================================================================================

ALLOCATE(MetricsAdv(3,0:PP_N,0:PP_N,0:PP_NZ,nElems,0:FV_ENABLED))
DO FVE=0,FV_ENABLED
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      MetricsAdv(1,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_fTilde(:,i,j,k,iElem,FVE))
      MetricsAdv(2,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_gTilde(:,i,j,k,iElem,FVE))
      MetricsAdv(3,i,j,k,iElem,FVE)=sJ(i,j,k,iElem,FVE)*NORM2(Metrics_hTilde(:,i,j,k,iElem,FVE))
    END DO; END DO; END DO
  END DO
END DO
END SUBROUTINE

!TODO compare to NS, adapt to 2D
!==================================================================================================================================
!> Calculate the time step for the current update of U for the Linear Scalar Advection Equation du/dt + a du/dx = 0
!==================================================================================================================================
FUNCTION CALCTIMESTEP(errType,doInit)
! MODULES
USE MOD_Globals
USE MOD_PreProc
#ifndef GNU
USE, INTRINSIC :: IEEE_ARITHMETIC,ONLY:IEEE_IS_NAN
#endif
USE MOD_Mesh_Vars    ,ONLY:sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,Elem_xGP,nElems
USE MOD_TimeDisc_Vars,ONLY:CFLScale,ViscousTimeStep,dtElem
USE MOD_Baseflow_Vars,ONLY:UBase
#if FV_ENABLED
USE MOD_FV_Vars
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL                         :: CalcTimeStep
INTEGER,INTENT(OUT)          :: errType
LOGICAL,INTENT(IN)           :: doInit
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: i,j,k,iElem
REAL,DIMENSION(PP_nVarBase)  :: UE
REAL                         :: TimeStepConv, TimeStepVisc, TimeStep(3)
REAL                         :: Max_Lambda(3),c,vsJ(3)
INTEGER                      :: FVE
!==================================================================================================================================
errType=0
CalcTimeStep=0.
IF(.NOT.doInit) RETURN

TimeStepConv=HUGE(1.)
TimeStepVisc=HUGE(1.)
DO iElem=1,nElems
  FVE = FV_Elems(iElem)
  Max_Lambda=0.
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    UE(:)=UBase(:,i,j,k,iElem)
    ! Convective Eigenvalues
    IF(IEEE_IS_NAN(UE(DENS)))THEN
      ERRWRITE(*,'(A,3ES16.7)')'Density NaN, Position= ',Elem_xGP(:,i,j,k,iElem)
      errType=1
    END IF
    c=UE(SOSP)
    vsJ=UE(VELV)*sJ(i,j,k,iElem,FVE)
    Max_Lambda(1)=MAX(Max_Lambda(1),ABS(SUM(Metrics_fTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(1,i,j,k,iElem,FVE))
    Max_Lambda(2)=MAX(Max_Lambda(2),ABS(SUM(Metrics_gTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(2,i,j,k,iElem,FVE))
    Max_Lambda(3)=MAX(Max_Lambda(3),ABS(SUM(Metrics_hTilde(:,i,j,k,iElem,FVE)*vsJ)) + &
                                              c*MetricsAdv(3,i,j,k,iElem,FVE))
  END DO; END DO; END DO ! i,j,k

  dtElem(iElem)=CFLScale(FVE)*2./SUM(Max_Lambda)
  TimeStepConv=MIN(TimeStepConv,dtElem(iElem))
  IF(IEEE_IS_NAN(TimeStepConv))THEN
    ERRWRITE(*,'(A,I0,A,I0)')'Convective timestep NaN on proc',myRank,' for element: ',iElem
    ERRWRITE(*,'(A,3ES16.7)')'Position: Elem_xGP(:1,1,1,iElem)=',Elem_xGP(:,1,1,1,iElem)
    ERRWRITE(*,*)'dt_conv=',TimeStepConv,' dt_visc=',TimeStepVisc
    errType=2
  END IF

END DO ! iElem=1,nElems

TimeStep(1)=TimeStepConv
TimeStep(2)=TimeStepVisc 
#if USE_MPI
TimeStep(3)=-errType ! reduce with timestep, minus due to MPI_MIN
CALL MPI_ALLREDUCE(MPI_IN_PLACE,TimeStep,3,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_FLEXI,iError)
errType=INT(-TimeStep(3))
#endif /*USE_MPI*/
ViscousTimeStep=(TimeStep(2) .LT. TimeStep(1))
CalcTimeStep=MINVAL(TimeStep(1:2))
END FUNCTION CALCTIMESTEP


!==================================================================================================================================
!> Timestep dummy routine. Since the LEE are linear, the timestep is only calculated during initialization.
!==================================================================================================================================
FUNCTION CALCTIMESTEP_dummy(errType)
USE MOD_Globals
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(OUT)          :: errType
!----------------------------------------------------------------------------------------------------------------------------------
REAL                         :: CalcTimeStep_dummy
!----------------------------------------------------------------------------------------------------------------------------------
errType=0
CalcTimeStep_dummy=TimeStep_global
END FUNCTION CALCTIMESTEP_dummy


!==================================================================================================================================
!> Deallocate CalcTimeStep arrays
!==================================================================================================================================
SUBROUTINE FinalizeCalctimestep()
! MODULES
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(MetricsAdv)
END SUBROUTINE

END MODULE MOD_CalcTimeStep
