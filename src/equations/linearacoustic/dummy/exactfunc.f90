!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
#include "eos.h"

!==================================================================================================================================
!> Routines providing initialization and initial solutions for the linear advection-diffusion equation
!==================================================================================================================================
MODULE MOD_Exactfunc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE DefineParametersExactFunc
  MODULE PROCEDURE DefineParametersExactFunc
END INTERFACE

INTERFACE InitExactFunc
  MODULE PROCEDURE InitExactFunc
END INTERFACE

INTERFACE ExactFunc
  MODULE PROCEDURE ExactFunc
END INTERFACE

INTERFACE CalcSource
  MODULE PROCEDURE CalcSource
END INTERFACE


PUBLIC::DefineParametersExactFunc
PUBLIC::InitExactFunc
PUBLIC::ExactFunc
PUBLIC::CalcSource
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersExactFunc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Exactfunc")
CALL prms%CreateIntFromStringOption(      'IniExactFunc', "Number of exactfunction to be used, to initialize the solution. "//&
                                                          "testcase,refstate,planewave,monopole,pulse2d,pulse3d,shearlayer")
CALL addStrListEntry('IniExactFunc','testcase' , -1)
CALL addStrListEntry('IniExactFunc','testcase' ,  0)
CALL addStrListEntry('IniExactFunc','refstate' ,  1)
CALL addStrListEntry('IniExactFunc','planewave',  2)
CALL addStrListEntry('IniExactFunc','monopole',   3)
CALL addStrListEntry('IniExactFunc','monopole3D', 31)
CALL addStrListEntry('IniExactFunc','pulse2d',    4)
CALL addStrListEntry('IniExactFunc','pulse3d',    5)
CALL prms%CreateRealArrayOption(    'xP',    "Pressure pulse initial position CASE(4,5) (x,y,z)")
CALL prms%CreateRealOption(         'sigma', "Initial pulse/source width  CASE(3,4,5)", '3')
CALL prms%CreateRealOption(         'frequency', "frequency for monopole source CASE(3)")
CALL prms%CreateRealOption(         'amplitude', "amplitude for monopole source CASE(3)")
END SUBROUTINE DefineParametersExactFunc

!==================================================================================================================================
!> Get some parameters needed for exact function
!==================================================================================================================================
SUBROUTINE InitExactFunc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools,   ONLY: GETINTFROMSTR,GETREAL,GETINT,GETREALARRAY
USE MOD_ExactFunc_Vars 
USE MOD_Equation_Vars, ONLY: IniExactFunc,IniRefState
! IMPLICIT VARIABLE HANDLING
 IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT EXACT FUNCTION...'

! Read in boundary parameters
IniExactFunc = GETINTFROMSTR('IniExactFunc')

! Read in parameters specific to certain init functions
SELECT CASE (IniExactFunc)
CASE DEFAULT
  IniRefState  = GETINT('IniRefState')
CASE(2) ! plane wave
CASE(3,31,4,5)
  IniRefState  = GETINT('IniRefState')
  xP = GETREALARRAY('xP',3,'(/0.,0.,0./)')
  sigmaSq = GETREAL('sigma','3')
  sigmaSq = sigmaSq**2
  SELECT CASE(IniExactFunc)
  CASE(3,31)
    f   = GETREAL('frequency','1')
    Amp = GETREAL('amplitude','1')
  END SELECT
END SELECT

SWRITE(UNIT_stdOut,'(A)')' INIT EXACT FUNCTION DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitExactFunc

!==================================================================================================================================
!> Specifies all the initial conditions. The state in conservative variables is returned.
!==================================================================================================================================
SUBROUTINE ExactFunc(ExactFunction,tIn,x,resu,RefStateOpt)
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars, ONLY: RefStateCons,IniRefState
USE MOD_Baseflow_Vars, ONLY: BaseState,varMeanFlow
USE MOD_Exactfunc_Vars,ONLY: xP,sigmaSq 
USE MOD_Timedisc_Vars, ONLY: fullBoundaryOrder,CurrentStage,dt,RKb,RKc,t
USE MOD_Testcase     , ONLY: ExactFuncTestcase
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: tIn                    !< input time (either time at RK stage or time at the beginning of
                                                          !< timestep if full boundary order is used (only with RK3)
REAL,INTENT(IN)                 :: x(3)                   !< coordinates to evaluate exact function
INTEGER,INTENT(IN)              :: ExactFunction          !< specifies the exact function to be used
REAL,INTENT(OUT)                :: Resu(PP_nVar)          !< output state in conservative variables
INTEGER,INTENT(IN),OPTIONAL     :: RefStateOpt            !< refstate to be used for exact func
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                            :: tEval
REAL                            :: Resu_t(PP_nVar),Resu_tt(PP_nVar) ! temporal state deriv in conservative variables
REAL                            :: Omega
REAL                            :: Cent(3),x0(3)
REAL                            :: Pi
REAL                            :: k,w
! case(4) Pulse
!==================================================================================================================================
tEval=MERGE(t,tIn,fullBoundaryOrder) ! prevent temporal order degradation, works only for RK3 time integration

Pi = ACOS(-1.)

Resu   =0.
Resu_t =0.
Resu_tt=0.

SELECT CASE (ExactFunction)
CASE DEFAULT
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)

CASE(0)
  CALL ExactFuncTestcase(tEval,x,Resu,Resu_t,Resu_tt)

CASE(1) ! constant
  Resu = RefStateCons(:,IniRefState)

CASE(2) ! 2D plane wave Birkefeld diss, p. 51
  k=2*pi
  w=(-1)*BaseState(6)*SQRT(k*k+k*k)
  Resu(1)=0.5*(BaseState(1)/BaseState(6))*SIN(k*x(1)+k*x(2)-w*tEval)
  Resu(2)=(-1)*0.5*(1/SQRT(2.0))*SIN(k*x(1)+k*x(2)-w*tEval)
  Resu(3)=(-1)*0.5*(1/SQRT(2.0))*SIN(k*x(1)+k*x(2)-w*tEval)
  Resu(4)=0.0
  Resu(5)=0.5*BaseState(1)*BaseState(6)*SIN(k*x(1)+k*x(2)-w*tEval)
  ! TODO: add Resu_t and Resu_tt for 2D plane wave case.
CASE(3,31) ! constant with monopole source term (Ewert & Schröder 2003, p.386)
  Resu = RefStateCons(:,IniRefState)
CASE(4)  ! 2D pulse
      Cent=x-xP
      Resu(5)=EXP(-LOG(2.)/sigmaSq*SUM(Cent(1:2)**2))
      Resu(1)=Resu(5)/BaseState(6)**2
CASE(5)  ! 3D pulse
      Cent=x-xP
      Resu(5)=EXP(-LOG(2.)/sigmaSq*SUM(Cent(1:3)**2))
      Resu(1)=Resu(5)/BaseState(6)**2
END SELECT ! ExactFunction

! For O3 LS 3-stage RK, we have to define proper time dependent BC
IF(fullBoundaryOrder)THEN ! add resu_t, resu_tt if time dependant
  SELECT CASE(CurrentStage)
  CASE(1)
    ! resu = g(t)
  CASE(2)
    ! resu = g(t) + dt/3*g'(t)
    Resu=Resu + dt*RKc(2)*Resu_t
  CASE(3,31)
    ! resu = g(t) + 3/4 dt g'(t) +5/16 dt^2 g''(t)
    Resu=Resu + RKc(3)*dt*Resu_t + RKc(2)*RKb(2)*dt*dt*Resu_tt
  CASE DEFAULT
    ! Stop, works only for 3 Stage O3 LS RK
    CALL abort(__STAMP__,&
               'Time-dependant exactfuntion works only for 3 Stage O3 LS RK!')
  END SELECT
END IF
END SUBROUTINE ExactFunc



!==================================================================================================================================
!> Compute source terms for some specific testcases and adds it to DG time derivative
!==================================================================================================================================
SUBROUTINE CalcSource(Ut,t)
! MODULES
USE MOD_Globals,       ONLY:Abort
USE MOD_Equation_Vars, ONLY:IniExactFunc,doCalcSource,readAcSourcesFromFile
USE MOD_Baseflow_Vars, ONLY:BaseState
USE MOD_ExactFunc_Vars,ONLY:f,xP,sigmaSq,Amp
USE MOD_Baseflow,      ONLY:CalcSourceBF
USE MOD_Equation_Vars, ONLY:useAcSources
USE MOD_AcSources,     ONLY:ApplyAcSources
USE MOD_Baseflow_Vars, ONLY:varMeanFlow
USE MOD_Mesh_Vars,     ONLY:nElems,Elem_xGP
USE MOD_PreProc
USE MOD_ApplyJacobianCons, ONLY:ApplyJacobianCons

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)     :: t                                       !< solution time
REAL,INTENT(INOUT)  :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems) !< solution time derivative
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k
REAL                :: Pi
REAL                :: source(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)
! case 3
REAL                :: Omega,sourceP,Cent(3)
!==================================================================================================================================
IF(varMeanFlow) THEN
  ! LEE source term arising from mean flow gradients...
  !CALL CalcSourceBF(source)
  source(:,:,:,:,:)=0.
ELSE
  source(:,:,:,:,:)=0.
END IF

IF(useAcSources .AND. readAcSourcesFromFile) CALL ApplyAcSources(source,t)

! Exact function dependent sources
Pi = ACOS(-1.)
SELECT CASE (IniExactFunc)
CASE(3) ! constant with 2D monopole source term (Ewert & Schröder 2003, p.386)
  Omega=2*pi*f  
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Cent=Elem_xGP(:,i,j,k,iElem)-xP
      sourceP=Amp*EXP(-LOG(2.)/sigmaSq*SUM(Cent(1:2)**2))*COS(Omega*t)
      source(5,i,j,k,iElem)=source(5,i,j,k,iElem)+sourceP
      source(1,i,j,k,iElem)=source(5,i,j,k,iElem)+sourceP/BaseState(6)**2
    END DO; END DO; END DO
  END DO
CASE(31) ! constant with 3D monopole source term
  Omega=2*pi*f  
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Cent=Elem_xGP(:,i,j,k,iElem)-xP
      sourceP=Amp*EXP(-LOG(2.)/sigmaSq*SUM(Cent(1:3)**2))*COS(Omega*t)
      source(5,i,j,k,iElem)=source(5,i,j,k,iElem)+sourceP
      source(1,i,j,k,iElem)=source(5,i,j,k,iElem)+sourceP/BaseState(6)**2
    END DO; END DO; END DO
  END DO
CASE DEFAULT
  ! No source -> do nothing and set marker to not run again
  IF((.NOT.varMeanFlow).AND.(.NOT.useAcSources)) THEN
    doCalcSource=.FALSE.
    RETURN
  END IF
END SELECT ! ExactFunction
CALL ApplyJacobianCons(source,toPhysical=.FALSE.)
Ut=Ut+source
END SUBROUTINE CalcSource

END MODULE MOD_Exactfunc
