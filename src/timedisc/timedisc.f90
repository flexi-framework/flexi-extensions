!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
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

!==================================================================================================================================
!> Module for the GTS Temporal discretization
!==================================================================================================================================
MODULE MOD_TimeDisc
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitTimeDisc
  MODULE PROCEDURE InitTimeDisc
END INTERFACE

INTERFACE TimeDisc
  MODULE PROCEDURE TimeDisc
END INTERFACE

INTERFACE FinalizeTimeDisc
  MODULE PROCEDURE FinalizeTimeDisc
END INTERFACE

PUBLIC :: InitTimeDisc,FinalizeTimeDisc
PUBLIC :: TimeDisc
PUBLIC :: DefineParametersTimeDisc
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersTimeDisc()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("TimeDisc")
CALL prms%CreateStringOption('TimeDiscMethod', "Specifies the type of time-discretization to be used, e.g. the name of&
                                               & a specific Runge-Kutta scheme. Possible values:\n"//&
                                               "  * standardrk3-3\n  * carpenterrk4-5\n  * niegemannrk4-14\n"//&
                                               "  * toulorgerk4-8c\n  * toulorgerk3-7c\n  * toulorgerk4-8f\n"//&
                                               "  * ketchesonrk4-20\n  * ketchesonrk4-18", value='CarpenterRK4-5')
CALL prms%CreateRealOption(  'TEnd',           "End time of the simulation (mandatory).")
CALL prms%CreateRealOption(  'CFLScale',       "Scaling factor for the theoretical CFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateRealOption(  'DFLScale',       "Scaling factor for the theoretical DFL number, typical range 0.1..1.0 (mandatory)")
CALL prms%CreateIntOption(   'maxIter',        "Stop simulation when specified number of timesteps has been performed.", value='-1')
CALL prms%CreateIntOption(   'NCalcTimeStepMax',"Compute dt at least after every Nth timestep.", value='1')
END SUBROUTINE DefineParametersTimeDisc

!==================================================================================================================================
!> Get information for end time and max time steps from ini file
!==================================================================================================================================
SUBROUTINE InitTimeDisc()
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars
USE MOD_ReadInTools         ,ONLY:GETREAL,GETINT,GETSTR
USE MOD_StringTools         ,ONLY:LowCase,StripSpaces
USE MOD_Overintegration_Vars,ONLY:NUnder
USE MOD_Filter_Vars         ,ONLY:NFilter,FilterType
USE MOD_Mesh_Vars           ,ONLY:nElems
USE MOD_IO_HDF5             ,ONLY:AddToElemData,ElementOut
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255):: TimeDiscMethod
INTEGER           :: NEff
!==================================================================================================================================
TimeDiscMethod = GETSTR('TimeDiscMethod','Carpenter RK4-5')
CALL StripSpaces(TimeDiscMethod)
CALL LowCase(TimeDiscMethod)

CALL SetTimeDiscCoefs(TimeDiscMethod)
SELECT CASE(TimeDiscType)
CASE('LSERKW2')
  TimeStep=>TimeStepByLSERKW2
CASE('LSERKK3')
  TimeStep=>TimeStepByLSERKK3
END SELECT

IF(TimeDiscInitIsDone)THEN
   SWRITE(*,*) "InitTimeDisc already called."
   RETURN
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT TIMEDISC...'

! Read the end time TEnd from ini file
TEnd     = GETREAL('TEnd')
! Read the normalized CFL number
CFLScale = GETREAL('CFLScale')
#if PARABOLIC
! Read the normalized DFL number
DFLScale = GETREAL('DFLScale')
#endif /*PARABOLIC*/
NEff=MIN(PP_N,NFilter,NUnder)
IF(FilterType.GT.2) NEff=PP_N!LAF,HESTHAVEN no timestep effect
CALL fillCFL_DFL(NEff,PP_N)
! Set timestep to a large number
dt=HUGE(1.)
! Read max number of iterations to perform
maxIter = GETINT('maxIter','-1')
nCalcTimeStepMax = GETINT('nCalcTimeStepMax','1')
SWRITE(UNIT_stdOut,'(A)') ' Method of time integration: '//TRIM(TimeDiscName)
ALLOCATE(dtElem(nElems))
dtElem=0.

CALL AddToElemData(ElementOut,'dt',dtElem)

TimediscInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT TIMEDISC DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitTimeDisc



!==================================================================================================================================
!> GTS Temporal discretization
!==================================================================================================================================
SUBROUTINE TimeDisc()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_TimeDisc_Vars       ,ONLY: TEnd,t,dt,tAnalyze,ViscousTimeStep,maxIter,Timestep,nRKStages,nCalcTimeStepMax,CurrentStage
USE MOD_Analyze_Vars        ,ONLY: Analyze_dt,WriteData_dt,tWriteData,nWriteData
USE MOD_AnalyzeEquation_Vars,ONLY: doCalcTimeAverage
USE MOD_Analyze             ,ONLY: Analyze
USE MOD_Equation_Vars       ,ONLY: StrVarNames
USE MOD_TestCase            ,ONLY: AnalyzeTestCase,CalcForcing
USE MOD_TestCase_Vars       ,ONLY: nAnalyzeTestCase,doTCSource
USE MOD_TimeAverage         ,ONLY: CalcTimeAverage
USE MOD_Restart_Vars        ,ONLY: DoRestart,RestartTime
USE MOD_CalcTimeStep        ,ONLY: CalcTimeStep
USE MOD_Output              ,ONLY: Visualize,PrintStatusLine
USE MOD_HDF5_Output         ,ONLY: WriteState,WriteBaseFlow
USE MOD_Mesh_Vars           ,ONLY: MeshFile,nGlobalElems
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U,JU
USE MOD_Overintegration     ,ONLY: Overintegration
USE MOD_Overintegration_Vars,ONLY: OverintegrationType
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_RecordPoints        ,ONLY: RecordPoints,WriteRP
USE MOD_RecordPoints_Vars   ,ONLY: RP_onProc
USE MOD_Sponge_Vars         ,ONLY: CalcPruettDamping
USE MOD_Indicator           ,ONLY: doCalcIndicator,CalcIndicator
#if FV_ENABLED
USE MOD_FV
#endif
#if GCL
USE MOD_GCL                 ,ONLY: GCLTimeDerivative_weakForm
#endif /*GCL*/
USE MOD_MoveMesh            ,ONLY: MoveMesh
use MOD_IO_HDF5
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                         :: dt_Min,dt_MinOld,dtAnalyze,dtEnd,tStart
INTEGER(KIND=8)              :: iter,iter_loc
REAL                         :: CalcTimeStart,CalcTimeEnd
INTEGER                      :: TimeArray(8)              !< Array for system time
INTEGER                      :: errType,nCalcTimestep,writeCounter
LOGICAL                      :: doAnalyze,doFinalize
!==================================================================================================================================

SWRITE(UNIT_StdOut,'(132("-"))')

! write number of grid cells and dofs only once per computation
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#GridCells : ',REAL(nGlobalElems)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs      : ',REAL(nGlobalElems)*REAL((PP_N+1)**PP_dim)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#Procs     : ',REAL(nProcessors)
SWRITE(UNIT_stdOut,'(A13,ES16.7)')'#DOFs/Proc : ',REAL(nGlobalElems*(PP_N+1)**PP_dim/nProcessors)


! Determine next write time, since it will be written into output file
t = MERGE(RestartTime,0.,DoRestart)
tWriteData=MIN(t+WriteData_dt,tEnd)
tAnalyze=MIN(t+Analyze_dt,tEnd)

! TODO: Should this be done before or after Overintegration? (see below)
! Write the state at time=0, i.e. the initial condition
!CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                      !FutureTime=tWriteData,isErrorFile=.FALSE.)

! --- Perform some preparational steps ---
! overintegrate solution first time
SELECT CASE(OverintegrationType)
CASE (1)
  CALL Overintegration(U)
CASE (2)
  CALL ApplyJacobianCons(U,toPhysical=.FALSE.,FVE=0)
  CALL Overintegration(U)
END SELECT

! Do first RK stage of first timestep to fill gradients
CurrentStage=1
CALL MoveMesh(t)
CALL ApplyJacobianCons(U,JU,toPhysical=.FALSE.) ! Transform initial solution to reference space
CALL DGTimeDerivative_weakForm(t)
#if GCL
CALL GCLTimeDerivative_weakForm()
#endif /*GCL*/
IF(doCalcIndicator) CALL CalcIndicator(U,t)

#if FV_ENABLED
! initial switch to FV sub-cells (must be called after DGTimeDerivative_weakForm, since indicator may require gradients)
IF(.NOT.DoRestart)THEN
  CALL FV_FillIni()
  CALL ApplyJacobianCons(U,JU,toPhysical=.FALSE.) ! Transform initial solution to reference space
END IF
#endif

IF(.NOT.DoRestart)THEN
  SWRITE(UNIT_StdOut,*)'WRITING INITIAL SOLUTION:'
ELSE
  SWRITE(UNIT_StdOut,*)'REWRITING SOLUTION:'
END IF

! TODO: Should this be done before or after Overintegration? (see above) For FV we need it after DGTimeDerivative_weakForm!
! Write the state at time=0, i.e. the initial condition
CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                      FutureTime=tWriteData,isErrorFile=.FALSE.)

CALL Visualize(t,U)

! No computation needed if tEnd=tStart!
IF((t.GE.tEnd).OR.maxIter.EQ.0) RETURN


tStart = t
iter=0
iter_loc=0
writeCounter=0
doAnalyze=.FALSE.
doFinalize=.FALSE.
! compute initial timestep
dt=CALCTIMESTEP(errType)
nCalcTimestep=0
dt_MinOld=-999.
IF(errType.NE.0) CALL abort(__STAMP__,&
#if EQNSYSNR == 3
  'Error: (1) density, (2) convective / (3) viscous timestep / muTilde (4) is NaN. Type/time:',errType,t)
#else
  'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
#endif

! Run initial analyze
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,*) 'Errors of initial solution:'
CALL Analyze(t,iter)
! fill recordpoints buffer (initialization/restart)
IF(RP_onProc) CALL RecordPoints(PP_nVar,StrVarNames,iter,t,.TRUE.)

CALL PrintStatusLine(t,dt,tStart,tEnd)

SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_StdOut,'(A,ES16.7)')'Initial Timestep  : ', dt
IF(ViscousTimeStep)THEN
  SWRITE(UNIT_StdOut,'(A)')' Viscous timestep dominates! '
END IF
#if FV_ENABLED
CALL FV_Info(1_8)
#endif
SWRITE(UNIT_StdOut,*)'CALCULATION RUNNING...'


! Run computation
CalcTimeStart=FLEXITIME()
DO
  CurrentStage=1
  IF(doCalcIndicator) CALL CalcIndicator(U,t)
#if FV_ENABLED
  CALL FV_Switch(U,AllowToDG=(nCalcTimestep.LT.1))
#endif
  CALL MoveMesh(t)                               !Calculate new mesh positions
  CALL ApplyJacobianCons(JU,U,toPhysical=.TRUE.) !Get current physical U 
  CALL DGTimeDerivative_weakForm(t)
#if GCL
  CALL GCLTimeDerivative_weakForm()
#endif /*GCL*/
  IF(nCalcTimestep.LT.1)THEN
    dt_Min=CALCTIMESTEP(errType)
    nCalcTimestep=MIN(FLOOR(ABS(LOG10(ABS(dt_MinOld/dt_Min-1.)**2.*100.+EPSILON(0.)))),nCalcTimeStepMax)
    dt_MinOld=dt_Min
    IF(errType.NE.0)THEN
      CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                            FutureTime=tWriteData,isErrorFile=.TRUE.)
      CALL abort(__STAMP__,&
     'Error: (1) density, (2) convective / (3) viscous timestep is NaN. Type/time:',errType,t)
    END IF
  END IF
  nCalcTimestep=nCalcTimestep-1

  dt=dt_Min
  dtAnalyze=HUGE(1.)
  IF(tAnalyze-t.LE.dt*(1.+1.E-4))THEN
    dtAnalyze=tAnalyze-t; doAnalyze=.TRUE.
  END IF
  dt=MIN(dt,dtAnalyze)

  dtEnd=HUGE(1.)
  IF(tEnd-t    .LE.dt*(1.+1.E-4))THEN
    dtEnd=tEnd-t;         doAnalyze=.TRUE.; doFinalize=.TRUE.
  END IF
  dt=MIN(dt,dtEnd)

  IF(doCalcTimeAverage) CALL CalcTimeAverage(.FALSE.,dt,t)
  IF(doTCSource) CALL CalcForcing(t,dt)


  CALL PrintStatusLine(t,dt,tStart,tEnd)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL TimeStep(t)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Perform Timestep using a global time stepping routine, attention: only RK3 has time dependent BC
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  iter=iter+1
  iter_loc=iter_loc+1
  t=t+dt
  IF(iter.EQ.maxIter)THEN
    tEnd=t; tAnalyze=t; tWriteData=t
    doAnalyze=.TRUE.; doFinalize=.TRUE.
  END IF

  ! Call DG operator to fill face data, fluxes, gradients for analyze
  IF(doAnalyze) THEN
    CALL MoveMesh(t)                               ! Move mesh to current position
    CALL ApplyJacobianCons(JU,U,toPhysical=.TRUE.) ! Get most recent physical state
    CALL DGTimeDerivative_weakForm(t)
  END IF

  ! Call your Analysis Routine for your Testcase here.
  IF((MOD(iter,nAnalyzeTestCase).EQ.0).OR.doAnalyze) CALL AnalyzeTestCase(t)
  ! evaluate recordpoints
  IF(RP_onProc) CALL RecordPoints(PP_nVar,StrVarNames,iter,t,doAnalyze)

  ! Analyze and output now
  IF(doAnalyze) THEN
    CalcTimeEnd=FLEXITIME()


    IF(MPIroot)THEN
      ! Get calculation time per DOF
      CalcTimeEnd=(CalcTimeEnd-CalcTimeStart)*REAL(nProcessors)/(REAL(nGlobalElems)*REAL((PP_N+1)**PP_dim)*REAL(iter_loc))/nRKStages
      CALL DATE_AND_TIME(values=TimeArray) ! get System time
      WRITE(UNIT_StdOut,'(132("-"))')
      WRITE(UNIT_stdOut,'(A,I2.2,A1,I2.2,A1,I4.4,A1,I2.2,A1,I2.2,A1,I2.2)') &
        ' Sys date   :    ',timeArray(3),'.',timeArray(2),'.',timeArray(1),' ',timeArray(5),':',timeArray(6),':',timeArray(7)
      WRITE(UNIT_stdOut,'(A,ES12.5,A)')' CALCULATION TIME PER STAGE/DOF: [',CalcTimeEnd,' sec ]'
      WRITE(UNIT_StdOut,'(A,ES16.7)')' Timestep   : ',dt_Min
      IF(ViscousTimeStep) WRITE(UNIT_StdOut,'(A)')' Viscous timestep dominates! '
      WRITE(UNIT_stdOut,'(A,ES16.7)')'#Timesteps  : ',REAL(iter)
    END IF !MPIroot
#if FV_ENABLED
    CALL FV_Info(iter_loc)
#endif

    ! Visualize data and write solution
    writeCounter=writeCounter+1
    IF((writeCounter.EQ.nWriteData).OR.doFinalize)THEN
      ! Write various derived data
      IF(doCalcTimeAverage) CALL CalcTimeAverage(.TRUE.,dt,t)
      IF(RP_onProc)         CALL WriteRP(PP_nVar,StrVarNames,t,.TRUE.)
      IF(CalcPruettDamping) CALL WriteBaseflow(TRIM(MeshFile),t)
      ! Write state file
      ! NOTE: this should be last in the series, so we know all previous data
      ! has been written correctly when the state file is present
      tWriteData=MIN(tAnalyze+WriteData_dt,tEnd)
      CALL WriteState(MeshFileName=TRIM(MeshFile),OutputTime=t,&
                            FutureTime=tWriteData,isErrorFile=.FALSE.)
      ! Visualize data
      CALL Visualize(t,U)
      writeCounter=0
    END IF

    ! do analysis
    CALL Analyze(t,iter)
    iter_loc=0
    CalcTimeStart=FLEXITIME()
    tAnalyze=  MIN(tAnalyze+Analyze_dt,  tEnd)
    doAnalyze=.FALSE.
  END IF

  IF(doFinalize) EXIT
END DO

END SUBROUTINE TimeDisc


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration: 2 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!> For moving meshes, we have to advance JU in time and calculate the new physical state after the time integration.
!> Naming of the variables is INCONSISTENT to get everything done with as little changes as possible.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKW2(t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U,JU,Ut,nTotalU
#if GCL
USE MOD_GCL                 ,ONLY: GCLTimeDerivative_weakForm
USE MOD_GCL_Vars            ,ONLY: nTotalGCL,Jac,Jac_t
USE MOD_Mesh_Vars           ,ONLY: sJ
#endif /*GCL*/
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_TimeDisc_Vars       ,ONLY: dt,RKA,RKb,RKc,nRKStages,CurrentStage
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_PruettDamping       ,ONLY: TempFilterTimeDeriv
USE MOD_Sponge_Vars         ,ONLY: CalcPruettDamping
USE MOD_Indicator           ,ONLY: doCalcIndicator,CalcIndicator
#if FV_ENABLED
USE MOD_FV                  ,ONLY: FV_Switch
USE MOD_FV_Vars             ,ONLY: FV_toDGinRK
#endif
USE MOD_MoveMesh            ,ONLY: MoveMesh
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: JUt_temp(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) ! temporal variable for JUt
#if GCL
REAL     :: Jact_temp(1,       0:PP_N,0:PP_N,0:PP_NZ,1:nElems) ! temporal variable for Jac_t
#endif
REAL     :: tStage,b_dt(1:nRKStages)
INTEGER  :: iStage
!===================================================================================================================================
! Premultiply with dt
b_dt=RKb*dt

IF(CalcPruettDamping) CALL TempFilterTimeDeriv(U,dt)

! First evaluation of DG operator already done in timedisc
CurrentStage=1
tStage=t
!CALL DGTimeDerivative_weakForm(tStage)      !allready called in timedisc
CALL VCopy(nTotalU,JUt_temp,Ut)                !JUt_temp = JUt
CALL VAXPBY(nTotalU,JU,Ut,ConstIn=b_dt(1))     !JU       = JU + JUt*b_dt(1)
#if GCL
#if FV_ENABLED
CALL CollectiveStop(__STAMP__, "FV and GCL are not validated. Think about the theory!!!")
#endif
!CALL GCLTimeDerivative_weakForm()
CALL VCopy(nTotalGCL,Jact_temp,Jac_t)               !Jact_temp = Jac_t
CALL VAXPBY(nTotalGCL,Jac,Jac_t,ConstIn=b_dt(1))    !Jac       = Jac + Jac_t*b_dt(1)
sJ(:,:,:,:,0) = 1./Jac(1,:,:,:,:)                   !Set current 1/J
#endif /*GCL*/

! Following steps
DO iStage=2,nRKStages
  CurrentStage=iStage
  tStage=t+dt*RKc(iStage)
  CALL MoveMesh(tStage)                          !Calculate new mesh positions
  CALL ApplyJacobianCons(JU,U,toPhysical=.TRUE.) !Get current physical U 
  IF(doCalcIndicator) CALL CalcIndicator(U,t)
#if FV_ENABLED
  CALL FV_Switch(U,JUt_temp,AllowToDG=FV_toDGinRK)
#endif
  CALL DGTimeDerivative_weakForm(tStage)
  CALL VAXPBY(nTotalU,JUt_temp,Ut,ConstOut=-RKA(iStage)) !JUt_temp = JUt - JUt_temp*RKA(iStage)
  CALL VAXPBY(nTotalU,JU,JUt_temp,ConstIn =b_dt(iStage)) !JU       = JU  + JUt_temp*b_dt(iStage)
#if GCL
  CALL GCLTimeDerivative_weakForm()
  CALL VAXPBY(nTotalGCL,Jact_temp,Jac_t,ConstOut=-RKA(iStage)) !Jact_temp = Jact - Jact_temp*RKA(iStage)
  CALL VAXPBY(nTotalGCL,Jac,Jact_temp,ConstIn =b_dt(iStage))   !Jac       = Jac + Jact_temp*b_dt(iStage)
  sJ(:,:,:,:,0) = 1./Jac(1,:,:,:,:)                            !Set current 1/J
#endif /*GCL*/
END DO
CurrentStage=1

END SUBROUTINE TimeStepByLSERKW2


!===================================================================================================================================
!> Low-Storage Runge-Kutta integration:  3 register version
!> This procedure takes the current time t, the time step dt and the solution at
!> the current time U(t) and returns the solution at the next time level.
!> RKA/b/c coefficients are low-storage coefficients, NOT the ones from butcher table.
!===================================================================================================================================
SUBROUTINE TimeStepByLSERKK3(t)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Vector
USE MOD_DG                  ,ONLY: DGTimeDerivative_weakForm
USE MOD_DG_Vars             ,ONLY: U,JU,Ut,nTotalU
#if GCL
USE MOD_GCL                 ,ONLY: GCLTimeDerivative_weakForm
USE MOD_GCL_Vars            ,ONLY: nTotalGCL,Jac,Jac_t
USE MOD_Mesh_Vars           ,ONLY: sJ
#endif /*GCL*/
USE MOD_ApplyJacobianCons   ,ONLY: ApplyJacobianCons
USE MOD_TimeDisc_Vars       ,ONLY: dt,RKdelta,RKg1,RKg2,RKg3,RKb,RKc,nRKStages,CurrentStage
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_PruettDamping       ,ONLY: TempFilterTimeDeriv
USE MOD_Sponge_Vars         ,ONLY: CalcPruettDamping
USE MOD_Indicator           ,ONLY: doCalcIndicator,CalcIndicator
#if FV_ENABLED
USE MOD_FV                  ,ONLY: FV_Switch
USE MOD_FV_Vars             ,ONLY: FV_toDGinRK
#endif
USE MOD_MoveMesh            ,ONLY: MoveMesh
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: t                                     !< current simulation time
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL     :: S2(    1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL     :: JUPrev(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#if GCL
REAL     :: Jac_S2( 1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
REAL     :: JacPrev(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)
#endif /*GCL*/
REAL     :: tStage,b_dt(1:nRKStages)
INTEGER  :: iStage
!===================================================================================================================================
IF(CalcPruettDamping) CALL TempFilterTimeDeriv(U,dt)


! Premultiply with dt
b_dt=RKb*dt

! Nomenclature:
! S1 == U, S2 == S2, S3 == UPrev

CurrentStage=1
tStage=t
CALL VCopy(nTotalU,JUprev,JU)                  !JUprev=JU
CALL VCopy(nTotalU,S2,JU)                      !S2=JU
!CALL DGTimeDerivative_weakForm(t)             ! allready called in timedisc
CALL DGTimeDerivative_weakForm(t)
CALL VAXPBY(nTotalU,JU,Ut,ConstIn=b_dt(1))     !JU      = JU + JUt*b_dt(1)
#if GCL
#if FV_ENABLED
CALL CollectiveStop(__STAMP__, "FV and GCL are not validated. Think about the theory!!!")
#endif
CALL VCopy(nTotalGCL,JacPrev,Jac)                !JacPrev=Jac
CALL VCopy(nTotalGCL,Jac_S2,Jac)                 !Jac_S2=Jac
!CALL GCLTimeDerivative_weakForm()
CALL VAXPBY(nTotalGCL,Jac,Jac_t,ConstIn=b_dt(1)) !Jac = Jac + Jac_t*b_dt(1)
sJ(:,:,:,:,0) = 1./Jac(1,:,:,:,:)                !Set current 1/J
#endif /*GCL*/

DO iStage=2,nRKStages
  CurrentStage=iStage
  tStage=t+dt*RKc(iStage)
  CALL MoveMesh(tStage)                          !Calculate new mesh positions
  CALL ApplyJacobianCons(JU,U,toPhysical=.TRUE.) !Get current physical U 
  IF(doCalcIndicator) CALL CalcIndicator(U,t)
#if FV_ENABLED
  CALL FV_Switch(U,JUprev,S2,AllowToDG=FV_toDGinRK)
#endif
  CALL DGTimeDerivative_weakForm(tStage)
  CALL VAXPBY(nTotalU,S2,JU,ConstIn=RKdelta(iStage))               !S2 = S2 + JU*RKdelta(iStage)
  CALL VAXPBY(nTotalU,JU,S2,ConstOut=RKg1(iStage),ConstIn=RKg2(iStage)) !JU = RKg1(iStage)*JU + RKg2(iStage)*S2
  CALL VAXPBY(nTotalU,JU,JUprev,ConstIn=RKg3(iStage))               !JU = JU + RKg3(ek)*JUprev
  CALL VAXPBY(nTotalU,JU,Ut,ConstIn=b_dt(iStage))                   !JU = JU + jUt*b_dt(iStage)
#if GCL
  CALL GCLTimeDerivative_weakForm()
  CALL VAXPBY(nTotalGCL,Jac_S2,Jac,ConstIn=RKdelta(iStage))                    !Jac_S2 = Jac_S2 + Jac*RKdelta(iStage)
  CALL VAXPBY(nTotalGCL,Jac,Jac_S2,ConstOut=RKg1(iStage),ConstIn=RKg2(iStage)) !Jac = RKg1(iStage)*Jac + RKg2(iStage)*Jac_S2
  CALL VAXPBY(nTotalGCL,Jac,JacPrev,ConstIn=RKg3(iStage))                      !Jac = Jac + RKg3(ek)*JacPrev
  CALL VAXPBY(nTotalGCL,Jac,Jac_t,ConstIn=b_dt(iStage))                        !Jac = Jac + Jac_t*b_dt(iStage)
  sJ(:,:,:,:,0) = 1./Jac(1,:,:,:,:)                                            !Set current 1/J
#endif /*GCL*/
END DO
CurrentStage=1

END SUBROUTINE TimeStepByLSERKK3


!===================================================================================================================================
!> Scaling of the CFL number, from paper GASSNER, KOPRIVA, "A comparision of the Gauss and Gauss-Lobatto
!> Discontinuous Galerkin Spectral Element Method for Wave Propagation Problems" .
!> For N=1-10 input CFLscale can now be (nearly) 1. and will be scaled adequately depending on
!> polynomial degree N, NodeType and TimeDisc method.
!===================================================================================================================================
SUBROUTINE FillCFL_DFL(Nin_CFL,Nin_DFL)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_TimeDisc_Vars,ONLY:CFLScale,CFLScaleAlpha
#if PARABOLIC
USE MOD_TimeDisc_Vars,ONLY:DFLScale,DFLScaleAlpha,RelativeDFL
#endif /*PARABOLIC*/
#if FV_ENABLED
USE MOD_TimeDisc_Vars,ONLY:CFLScaleFV
#endif /*FV*/
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: Nin_CFL !< input polynomial degree for advection terms
INTEGER,INTENT(IN) :: Nin_DFL !< input polynomial degree for viscous terms
                              !< for overintegration via Filtering, the gradients and thus the visscous flux remains of order N+1
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: alpha
#if !(PARABOLIC)
INTEGER            :: dummy
#endif
!===================================================================================================================================
! CFL in DG depends on the polynomial degree
! Runge-Kutta methods
alpha    = CFLScaleAlpha(MIN(15,Nin_CFL))
CFLScale(0) = CFLScale(0)*alpha
#if FV_ENABLED
CFLScale(1) = CFLScale(1)*CFLScaleFV/(PP_N+1.) ! equidistant distribution
#endif
IF((Nin_CFL.GT.15).OR.(CFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_StdOut,'(132("!"))')
  SWRITE(UNIT_StdOut,'(A)')'Warning: The chosen CFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_StdOut,'(132("!"))')
END IF
!scale with 2N+1
CFLScale(0) = CFLScale(0)/(2.*Nin_CFL+1.)
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   CFL (DG/FV):',CFLScale

#if PARABOLIC
!########################### DFL ########################################
! DFL in DG depends on the polynomial degree
! since DFl is only on real axis, stability numbers are defined for RK3 and then scaled for RK4

alpha = DFLScaleAlpha(MIN(10,Nin_DFL))*RelativeDFL
DFLScale(0)=DFLScale(0)*alpha
#if FV_ENABLED
DFLScale(1)=DFLScale(1)*DFLScaleAlpha(1)*RelativeDFL/(PP_N+1.)**2
#endif
IF((Nin_DFL.GT.10).OR.(DFLScale(0).GT.alpha))THEN
  SWRITE(UNIT_StdOut,'(132("!"))')
  SWRITE(UNIT_StdOut,'(A)')'Warning: The chosen DFL number may be too high for the selected polynomial degree!'
  SWRITE(UNIT_StdOut,'(132("!"))')
END IF
DFLScale(0) = DFLScale(0)/(2.*Nin_DFL+1.)**2
SWRITE(UNIT_stdOut,'(A,2ES16.7)') '   DFL (DG/FV):',DFLScale
#else
dummy = Nin_DFL ! prevent compile warning
#endif /*PARABOLIC*/
END SUBROUTINE fillCFL_DFL

!==================================================================================================================================
!> Finalizes variables necessary for timedisc subroutines
!==================================================================================================================================
SUBROUTINE FinalizeTimeDisc()
! MODULES
USE MOD_TimeDisc_Vars
IMPLICIT NONE
!==================================================================================================================================
TimeDiscInitIsDone = .FALSE.
SDEALLOCATE(dtElem)
SDEALLOCATE(RKA)
SDEALLOCATE(RKb)
SDEALLOCATE(RKc)
SDEALLOCATE(RKg1)
SDEALLOCATE(RKg2)
SDEALLOCATE(RKg3)
SDEALLOCATE(RKdelta)
NULLIFY(TimeStep)
END SUBROUTINE FinalizeTimeDisc


END MODULE MOD_TimeDisc
