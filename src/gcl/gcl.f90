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

!===================================================================================================================================
!> \brief Routines concerning the Geometric Conservation Law (GCL).
!> 
!> Contains the initialization and finalization of the GCL global variables and the routine to 
!> compute the GCL spatial operator/residual (Jac_t) from the mesh velocity.
!===================================================================================================================================
MODULE MOD_GCL
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------

! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE InitGCL
  MODULE PROCEDURE InitGCL
END INTERFACE

INTERFACE GCLTimeDerivative_weakForm
  MODULE PROCEDURE GCLTimeDerivative_weakForm
END INTERFACE

INTERFACE FinalizeGCL
  MODULE PROCEDURE FinalizeGCL
END INTERFACE

#ifdef SPLIT_DG
PUBLIC::DefineParametersGCL
#endif /*SPLIT_DG*/
PUBLIC::InitGCL,GCLTimeDerivative_weakForm,FinalizeGCL
!===================================================================================================================================

CONTAINS

#ifdef SPLIT_DG
!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersGCL()
! MODULES
USE MOD_ReadInTools,        ONLY: prms
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("GCL")
CALL prms%CreateLogicalOption('GCLSplitForm',"Switch if the GCL should be discretized using a split form or not.",'.TRUE.')
END SUBROUTINE DefineParametersGCL
#endif /*SPLIT_DG*/

!===================================================================================================================================
!> Allocate global variable Jac (Jacobian Determinant, conserved quantity of the GCL) and Jac_t (GCL time derivative).
!> Also GCLFlux (Fluxes over sides) are allocated.
!> All arrays have a first dimension of 1 so normal FLEXI routines can be used.
!===================================================================================================================================
SUBROUTINE InitGCL()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_GCL_Vars
USE MOD_ReadInTools,        ONLY: GETLOGICAL
USE MOD_GCL_VolInt,         ONLY: GCL_VolInt_pointer,GCL_VolInt_weakForm
#ifdef SPLIT_DG
USE MOD_GCL_VolInt,         ONLY: GCL_VolInt_splitForm
#endif /*SPLIT_DG*/
USE MOD_DG_Vars,            ONLY: nDOFElem
USE MOD_Restart_Vars,       ONLY: RestartInitIsDone
USE MOD_Interpolation_Vars, ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars,          ONLY: MeshInitIsDone,nSides,sJ,nElems
USE MOD_IO_HDF5,            ONLY: AddToFieldData,FieldOut
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=8)                :: DataSetName
CHARACTER(LEN=255)              :: VarNames(1)
INTEGER                         :: nVal(4)
!===================================================================================================================================
IF((.NOT.InterpolationInitIsDone).OR.(.NOT.MeshInitIsDone).OR.(.NOT.RestartInitIsDone).OR.GCLInitIsDone)THEN
   CALL CollectiveStop(__STAMP__,'InitGCL not ready to be called or already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT GCL...'

! A procedure pointer is used for the volume integral routine, so we can switch between the split version
! of the GCL and the standard formulation if the code is compiled with the split option.
#ifdef SPLIT_DG
IF (GETLOGICAL('GCLSplitForm','.TRUE.')) THEN
  GCL_Volint_pointer=>GCL_VolInt_splitForm
ELSE
  GCL_Volint_pointer=>GCL_VolInt_weakForm
END IF
#else
GCL_Volint_pointer=>GCL_VolInt_weakForm
#endif /*SPLIT_DG*/

! the local Jacobian Determinant solution
ALLOCATE(Jac(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
! set the inital value, computed during Mesh-Init
Jac(1,:,:,:,:) = 1./sJ(:,:,:,:,0)
! the GCL time derivative computed with the DG scheme
ALLOCATE(Jac_t(1,0:PP_N,0:PP_N,0:PP_NZ,nElems))
Jac_t = 0.

! unique flux per side
ALLOCATE(GCLFlux_master(1,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(GCLFlux_slave( 1,0:PP_N,0:PP_NZ,1:nSides))
GCLFlux_master=0.
GCLFlux_slave =0.

! variables for performance tricks
nTotalGCL=1*nDOFElem*nElems

! Add Jacobian computed by GCL to HDF5 output
DataSetName = 'Jacobian'
VarNames(1) = 'Jacobian_GCL'
nVal        = (/1,PP_N+1,PP_N+1,PP_NZ+1/)
CALL AddToFieldData(FieldOut,nVal,DataSetName,VarNames,Jac,doSeparateOutput=.TRUE.)

GCLInitIsDone=.TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT GCL DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitGCL


!===================================================================================================================================
!> Computes the GCL time derivative consisting of volume integral and surface integral for the whole field.
!===================================================================================================================================
SUBROUTINE GCLTimeDerivative_weakForm()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Vector
USE MOD_GCL_Vars,      ONLY: Jac_t,GCLFlux_master,GCLFlux_slave,nTotalGCL
USE MOD_GCL_SurfInt,   ONLY: GCL_SurfInt
USE MOD_GCL_VolInt,    ONLY: GCL_VolInt_pointer
USE MOD_GCL_FillFlux,  ONLY: GCL_FillFlux
USE MOD_DG_Vars,       ONLY: L_HatPlus,L_HatMinus
USE MOD_FillMortarMesh,ONLY: Flux_MortarGCL
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,           ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,     ONLY: nSides
#endif
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
! Nullify the time derivative of the GCL
CALL VNullify(nTotalGCL,Jac_t)

#if USE_MPI
! fill MPI side surface fluxes and send them  to the slave procs
CALL StartReceiveMPIData(GCLFlux_slave,DataSizeSideGCL,1,nSides,MPIRequest_Flux(:,SEND),SendID=1) ! Receive MINE
CALL GCL_FillFlux(GCLFlux_master,GCLFlux_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(GCLFlux_slave,DataSizeSideGCL,1,nSides,MPIRequest_Flux(:,RECV),SendID=1) ! Send YOUR
#endif /* MPI*/

! compute volume integral contribution and add to Jac_t
CALL GCL_VolInt_pointer(Jac_t)

! fill the interior surface fluxes
CALL GCL_FillFlux(GCLFlux_master,GCLFlux_slave,doMPISides=.FALSE.)
! Project fluxes from small mortar sides to big mortar sides
CALL Flux_MortarGCL(GCLFlux_master,GCLFlux_slave,doMPISides=.FALSE.,weak=.TRUE.)
! compute surface integral contribution and add to Jac_t
CALL GCL_SurfInt(PP_N,GCLFlux_master,GCLFlux_slave,Jac_t,.FALSE.,L_HatMinus,L_HatPlus)

#if USE_MPI
! Complete send / receive
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_Flux)
! Project fluxes from small mortar sides to big mortar sides
CALL Flux_MortarGCL(GCLFlux_master,GCLFlux_slave,doMPISides=.TRUE.,weak=.TRUE.)
! Surface integral for the MPI sides
CALL GCL_SurfInt(PP_N,GCLFlux_master,GCLFlux_slave,Jac_t,.TRUE.,L_HatMinus,L_HatPlus)
#endif

! Bring time derivative on other side
Jac_t = -Jac_t

END SUBROUTINE GCLTimeDerivative_weakForm

!===================================================================================================================================
!> Deallocate global variable Jac (Jacobian Determinant) ,Jac_t and GCL Fluxes.
!===================================================================================================================================
SUBROUTINE FinalizeGCL()
! MODULES
USE MOD_GCL_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!===================================================================================================================================
SDEALLOCATE(Jac)
SDEALLOCATE(Jac_t)
SDEALLOCATE(GCLFlux_master)
SDEALLOCATE(GCLFlux_slave)
GCLInitIsDone = .FALSE.
END SUBROUTINE FinalizeGCL

END MODULE MOD_GCL
