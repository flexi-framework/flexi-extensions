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
!> Routines providing initialization and initial solutions for the linearized euler equation
!==================================================================================================================================
MODULE MOD_Equation
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitEquation
  MODULE PROCEDURE InitEquation
END INTERFACE

INTERFACE GetPrimitiveStateSurface
  MODULE PROCEDURE GetPrimitiveStateSurface
END INTERFACE

INTERFACE GetConservativeStateSurface
  MODULE PROCEDURE GetConservativeStateSurface
END INTERFACE

INTERFACE FinalizeEquation
  MODULE PROCEDURE FinalizeEquation
END INTERFACE


PUBLIC:: DefineParametersEquation,InitEquation,FinalizeEquation
PUBLIC:: GetPrimitiveStateSurface,GetConservativeStateSurface
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersEquation()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_Riemann     ,ONLY: DefineParametersRiemann
USE MOD_Baseflow    ,ONLY: DefineParametersBaseFlow
USE MOD_Flux        ,ONLY: DefineParametersFlux
USE MOD_AcSources   ,ONLY: DefineParametersAcSources,DefineParametersAcMask
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("Equation")
CALL prms%CreateIntOption(      'IniRefState',  "Refstate required for initialization.")
CALL prms%CreateRealArrayOption('RefState',     "State(s) in primitive variables (density, velx, vely, velz, pressure).",&
                                                multiple=.TRUE.)
CALL prms%CreateRealOption(     'kappa',        "Heat capacity ratio / isentropic exponent", '1.4')
!useAcSources:
!when true swith for:
!1. InitEquation: InitAcSources()
!2. timedisc: prepAcSources(t,dt)
!3. DGTimeDerivative_weakForm => CalcSource(Ut,t) => ApplyAcSources(source,t) 
CALL prms%CreateLogicalOption(  'useAcSources', "Addition of read-in acoustic sources", 'F')
CALL prms%CreateLogicalOption(  'readAcSourcesFromFile', "Read the acoustic sources from files", 'T')
CALL DefineParametersRiemann()
CALL DefineParametersBaseFlow()
CALL DefineParametersFlux()
CALL DefineParametersAcSources()
CALL DefineParametersAcMask()

END SUBROUTINE DefineParametersEquation

!==================================================================================================================================
!> Read equation parameters (advection velocity, diffusion coeff, exact function)  from the ini file
!==================================================================================================================================
SUBROUTINE InitEquation()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars
USE MOD_Exactfunc         ,ONLY: InitExactFunc
USE MOD_Baseflow          ,ONLY: InitBaseflow
USE MOD_AcSources         ,ONLY: InitAcSources
USE MOD_Riemann           ,ONLY: InitRiemann
USE MOD_GetBoundaryFlux,   ONLY: InitBC
USE MOD_Testcase          ,ONLY: InitTestcase
USE MOD_ReadInTools       ,ONLY: CountOption,GETREALARRAY,GETSTR,GETREAL,GETLOGICAL
USE MOD_CalcTimeStep      ,ONLY: InitCalctimestep
USE MOD_Mesh_Vars         ,ONLY: nElems,nSides,Elem_xGP,MeshInitIsDone
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Flux              ,ONLY: InitFlux
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.EquationInitIsDone.OR.(.NOT.MeshInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitLEE not ready to be called or already called.")
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT LINEARIZED EULER EQUATIONS...'

s43=4./3.
s23=2./3.

! Always set docalcsource true, set false by calcsource itself on first run if not needed
doCalcSource=.TRUE.

! Gas constants
Kappa    =GETREAL('kappa')
KappaM1  =Kappa-1.
sKappaM1 =1./KappaM1
KappaP1  =Kappa+1.
sKappaP1 =1./(KappaP1)

IniRefState  = 0

! Read Boundary information / RefStates / perform sanity check
nRefState=CountOption('RefState')
IF(IniRefState.GT.nRefState)THEN
  CALL CollectiveStop(__STAMP__,&
    'ERROR: Ini not defined! (Ini,nRefState):',IniRefState,REAL(nRefState))
END IF

IF(nRefState .GT. 0)THEN
  ALLOCATE(RefStateCons(PP_nVar,nRefState))
  DO i=1,nRefState
    RefStateCons(1:5,i)  = GETREALARRAY('RefState',5)
  END DO
END IF

ALLOCATE(source(1:PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems))

! Initialize Base flow
CALL InitBaseflow()
source=0.

! Initialize acoustic source terms
useAcSources=GETLOGICAL('useAcSources','F')
readAcSourcesFromFile=GETLOGICAL('readAcSourcesFromFile','F')
IF(useAcSources .OR. readAcSourcesFromFile) CALL InitAcSources()

! Call initialization of exactfunc
CALL InitExactFunc()

! Initialize flux 
CALL InitFlux()

! Initialize Riemann solvers to be in volume and on BCs
CALL InitRiemann()

! Initialize current testcase
CALL InitTestcase()

! Initialize timestep calculation
CALL InitCalctimestep()

CALL InitBC()

EquationInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT LEE DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitEquation


!==================================================================================================================================
!> Converts conservative solution vector to primitive variables
!==================================================================================================================================
SUBROUTINE GetPrimitiveStateSurface(U_master,U_slave,UPrim_master,UPrim_slave)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY:nSides
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_master(        PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(IN)  :: U_slave(         PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
REAL,INTENT(OUT) :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(OUT) :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Copy the coservative state to the primitive arrays
UPrim_slave = U_slave
UPrim_master = U_master
END SUBROUTINE GetPrimitiveStateSurface


!==================================================================================================================================
!> Converts primite solution vector to conservative variables
!==================================================================================================================================
SUBROUTINE GetConservativeStateSurface(UPrim_master,UPrim_slave,U_master,U_slave, mask_master, mask_slave, mask_ref)
! MODULES
USE MOD_Preproc
USE MOD_Mesh_Vars,ONLY: firstInnerSide,firstMPISide_YOUR,lastMPISide_YOUR,nSides
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on master sides
REAL,INTENT(IN)    :: UPrim_slave( PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< primitive solution on slave sides
REAL,INTENT(OUT)   :: U_master(        PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on master sides
REAL,INTENT(OUT)   :: U_slave(         PP_nVar,0:PP_N,0:PP_NZ,1:nSides) !< conservative solution on slave sides
INTEGER,INTENT(IN) :: mask_master(1:nSides)                            !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_slave (1:nSides)                            !< mask: only convert solution if mask(SideID) == mask_ref
INTEGER,INTENT(IN) :: mask_ref                                         !< reference value for mask comparison
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Copy the coservative state to the primitive arrays
U_master = UPrim_master
U_slave  = UPrim_slave 
END SUBROUTINE GetConservativeStateSurface

!==================================================================================================================================
!> Finalizes the equation
!==================================================================================================================================
SUBROUTINE FinalizeEquation()
! MODULES
USE MOD_Equation_Vars
USE MOD_Testcase        ,ONLY: FinalizeTestcase
USE MOD_Riemann         ,ONLY: FinalizeRiemann
USE MOD_CalcTimeStep    ,ONLY: FinalizeCalctimestep
USE MOD_Baseflow        ,ONLY: FinalizeBaseflow
USE MOD_GetBoundaryFlux, ONLY: FinalizeBC
USE MOD_AcSources       ,ONLY: FinalizeAcSources
IMPLICIT NONE
!==================================================================================================================================
CALL FinalizeTestcase()
CALL FinalizeRiemann()
CALL FinalizeCalctimestep()
CALL FinalizeBaseflow()
!TODO vgl. Navier Stokes CALL FinalizeBC()
CALL FinalizeAcSources()
CALL FinalizeBC()
SDEALLOCATE(source)
EquationInitIsDone = .FALSE.
END SUBROUTINE FinalizeEquation

END MODULE MOD_Equation
