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
!> Routines to provide boundary conditions for the domain. Fills the boundary part of the fluxes list.
!==================================================================================================================================
MODULE MOD_GetBoundaryFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part --------------------------------------------------------------------------------------------------------------------
INTERFACE InitBC
  MODULE PROCEDURE InitBC
END INTERFACE

INTERFACE GetBoundaryFlux
  MODULE PROCEDURE GetBoundaryFlux
END INTERFACE

#if FV_ENABLED
INTERFACE GetBoundaryFVgradient
  MODULE PROCEDURE GetBoundaryFVgradient
END INTERFACE
#endif

INTERFACE FinalizeBC
  MODULE PROCEDURE FinalizeBC
END INTERFACE

PUBLIC :: InitBC
PUBLIC :: GetBoundaryFlux
#if FV_ENABLED
PUBLIC :: GetBoundaryFVgradient
#endif
PUBLIC :: FinalizeBC
!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Initialize boundary conditions. Read parameters and sort boundary conditions by types.
!> Call boundary condition specific init routines.
!==================================================================================================================================
SUBROUTINE InitBC()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars     ,ONLY: EquationInitIsDone
USE MOD_Equation_Vars     ,ONLY: BCData,nBCByType,BCSideID
USE MOD_Interpolation_Vars,ONLY: InterpolationInitIsDone
USE MOD_Mesh_Vars         ,ONLY: MeshInitIsDone,nBCSides,BC,BoundaryType,nBCs
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i,iSide
INTEGER :: locType,locState
INTEGER :: MaxBCState,MaxBCStateGlobal
!==================================================================================================================================
IF((.NOT.InterpolationInitIsDone).AND.(.NOT.MeshInitIsDone).AND.(.NOT.EquationInitIsDone))THEN
  CALL CollectiveStop(__STAMP__,&
    "InitBC not ready to be called or already called.")
END IF
! determine globally max MaxBCState
MaxBCState = 0
DO iSide=1,nBCSides
  locType =BoundaryType(BC(iSide),BC_TYPE)
  locState=BoundaryType(BC(iSide),BC_STATE)
END DO
MaxBCStateGLobal=MaxBCState
#if USE_MPI
CALL MPI_ALLREDUCE(MPI_IN_PLACE,MaxBCStateGlobal,1,MPI_INTEGER,MPI_MAX,MPI_COMM_FLEXI,iError)
#endif /*USE_MPI*/

! Allocate buffer array to store temp data for all BC sides
ALLOCATE(BCData(PP_nVar,0:PP_N,0:PP_NZ,nBCSides))
BCData=0.

! Initialize State File Boundary condition
DO i=1,nBCs
  locType =BoundaryType(i,BC_TYPE)
END DO

! Count number of sides of each boundary
ALLOCATE(nBCByType(nBCs))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i) nBCByType(i)=nBCByType(i)+1
  END DO
END DO

! Sort BCs by type, store SideIDs
ALLOCATE(BCSideID(nBCs,MAXVAL(nBCByType)))
nBCByType=0
DO iSide=1,nBCSides
  DO i=1,nBCs
    IF(BC(iSide).EQ.i)THEN
      nBCByType(i)=nBCByType(i)+1
      BCSideID(i,nBCByType(i))=iSide
    END IF
  END DO
END DO

END SUBROUTINE InitBC


!==================================================================================================================================
!> Computes the boundary values for a given Cartesian mesh face (defined by FaceID)
!> BCType: 1...periodic, 2...exact BC
!==================================================================================================================================
SUBROUTINE GetBoundaryFlux(SideID,t,Nloc,Flux,U_master,UB,                &
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
USE MOD_PreProc
USE MOD_Globals      ,ONLY: Abort
USE MOD_Mesh_Vars    ,ONLY: BC,nBCs,BoundaryType,nSides
USE MOD_Exactfunc    ,ONLY: ExactFunc
USE MOD_Equation_Vars,ONLY: IniExactFunc
USE MOD_Equation_Vars,ONLY: nBCByType,BCSideID
USE MOD_Flux         ,ONLY: EvalEulerFlux1D
USE MOD_Riemann      ,ONLY: Riemann
#if FV_ENABLED
USE MOD_FV_Vars      ,ONLY: FV_Elems_master
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)                   :: SideID  
REAL,INTENT(IN)                      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)                   :: Nloc    !< polynomial degree
REAL,INTENT(IN)                      :: U_master( PP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides) !< inner surface solution
REAL,INTENT(IN)                      :: UB( PP_nVarBase,0:Nloc,0:ZDIM(Nloc),1:nSides) !< baseflow state
REAL,DIMENSION(      3,0:Nloc,0:ZDIM(Nloc),0:FV_ENABLED,1:nSides),INTENT(IN)  :: NormVec,TangVec1,TangVec2
REAL,DIMENSION(      3,0:Nloc,0:ZDIM(Nloc),0:FV_ENABLED,1:nSides),INTENT(IN)  :: Face_xGP   !< positions of surface flux points
REAL,DIMENSION(PP_nVar,0:Nloc,0:ZDIM(Nloc),1:nSides),INTENT(OUT) :: Flux       !< resulting boundary fluxes
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                              :: iBC,iSide,p,q
INTEGER                              :: BCType,BCState,nBCLoc
REAL,DIMENSION(PP_nVar,0:Nloc,0:Nloc):: U_Face_loc
REAL,DIMENSION(PP_nVar)              :: U_loc,F_loc
REAL,DIMENSION(PP_nVarBase)          :: UB_loc
INTEGER                              :: FVEM=0 ! FV_Elems_master
!==================================================================================================================================
DO iBC=1,nBCs
  !IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(BC(SideID),BC_TYPE)
  BCState=BoundaryType(BC(SideID),BC_STATE)
  !nBCLoc =nBCByType(iBC)

  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) !Exact function or refstate
    ! BCState specifies refstate to be used, if 0 then use iniexactfunc
    !DO iSide=1,nBCLoc
#if FV_ENABLED
      FVEM = FV_Elems_master(SideID)
#endif
      !IF(BCState.EQ.0)THEN
        DO q=0,ZDIM(Nloc); DO p=0,Nloc
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,FVEM,SideID),U_Face_loc(:,p,q))
        END DO; END DO
      !ELSE
        !DO q=0,Nloc; DO p=0,Nloc
          !U_Face_loc(:,p,q) = RefStateCons(:,BCState)
        !END DO; END DO
      !END IF
    CALL Riemann(Nloc,Flux(:,:,:,SideID),U_master(:,:,:,SideID),U_Face_loc,UB(:,:,:,SideID),UB(:,:,:,SideID),&
                 NormVec(:,:,:,FVEM,SideID),TangVec1(:,:,:,FVEM,SideID),TangVec2(:,:,:,FVEM,SideID),.TRUE.)

  CASE(3) ! Wall BC, variant mirrored normal velocity
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      
      ! transform state to normal system to get mirror state
      U_loc(1)=U_master(1,p,q,SideID) 
      U_loc(5)=U_master(5,p,q,SideID) 
      U_loc(2)=(-1)*DOT_PRODUCT(U_master(2:4,p,q,SideID),NormVec(:,p,q,FVEM,SideID))
      U_loc(3)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      U_loc(4)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! transform baseflow to normal system
      UB_loc(:)=UB(:,p,q,SideID) 
      UB_loc(2)=0.
      UB_loc(3)=DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      UB_loc(4)=DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! get local flux of mirror state
      CALL EvalEulerFlux1D(U_loc,UB_loc,F_loc)
      ! rotate flux back to global coordinate system
      Flux(DENS,p,q,SideID)=F_loc(DENS)
      Flux(VELV,p,q,SideID)= NormVec(:,p,q,FVEM,SideID)*F_loc(VEL1)  &
                            +TangVec1(:,p,q,FVEM,SideID)*F_loc(VEL2) + TangVec2(:,p,q,FVEM,SideID)*F_loc(VEL3)
      Flux(PRES,p,q,SideID)=F_loc(PRES)
    END DO; END DO
  CASE(4) ! LEE wall BC, variant zero normal velocity
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! transform state to normal system and set normal component to zero
      U_loc(1)=U_master(1,p,q,SideID) 
      U_loc(5)=U_master(5,p,q,SideID) 
      U_loc(2)=0.
      U_loc(3)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      U_loc(4)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! transform baseflow to normal system
      UB_loc(:)=UB(:,p,q,SideID) 
      UB_loc(2)=0.
      UB_loc(3)=DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      UB_loc(4)=DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! get local flux of mirror state
      CALL EvalEulerFlux1D(U_loc,UB_loc,F_loc)
      ! rotate flux back to global coordinate system
      Flux(DENS,p,q,SideID)=F_loc(DENS)
      Flux(VELV,p,q,SideID)= NormVec(:,p,q,FVEM,SideID)*F_loc(VEL1)  &
                            +TangVec1(:,p,q,FVEM,SideID)*F_loc(VEL2) + TangVec2(:,p,q,FVEM,SideID)*F_loc(VEL3)
      Flux(PRES,p,q,SideID)=F_loc(PRES)
    END DO; END DO
  CASE(41) ! LEE wall BC, variant zero normal velocity
    DO q=0,ZDIM(Nloc); DO p=0,Nloc
      ! transform state to normal system and set normal component to zero
      U_loc(1)=U_master(1,p,q,SideID) 
      U_loc(5)=U_master(5,p,q,SideID) 
      U_loc(2)=0.
      U_loc(3)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      U_loc(4)=DOT_PRODUCT(U_master(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! transform baseflow to normal system
      UB_loc(:)=UB(:,p,q,SideID) 
      UB_loc(2)=0.
      UB_loc(3)=0. ! DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec1(:,p,q,FVEM,SideID))
      UB_loc(4)=0. ! DOT_PRODUCT(UB(2:4,p,q,SideID),TangVec2(:,p,q,FVEM,SideID))

      ! get local flux of mirror state
      CALL EvalEulerFlux1D(U_loc,UB_loc,F_loc)
      ! rotate flux back to global coordinate system
      Flux(DENS,p,q,SideID)=F_loc(DENS)
      Flux(VELV,p,q,SideID)= NormVec(:,p,q,FVEM,SideID)*F_loc(VEL1)  &
                            +TangVec1(:,p,q,FVEM,SideID)*F_loc(VEL2) + TangVec2(:,p,q,FVEM,SideID)*F_loc(VEL3)
      Flux(PRES,p,q,SideID)=F_loc(PRES)
    END DO; END DO

     
  
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryflux.f90!')
  END SELECT ! BCType
END DO

END SUBROUTINE GetBoundaryFlux

#if FV_ENABLED
SUBROUTINE GetBoundaryFVgradient(t,gradU,U_inner)
!==================================================================================================================================
!==================================================================================================================================
! MODULES
USE MOD_Globals       ,ONLY: Abort
USE MOD_PreProc
USE MOD_Exactfunc     ,ONLY: ExactFunc
USE MOD_Equation_Vars ,ONLY: IniExactFunc
USE MOD_Equation_Vars ,ONLY: nBCByType,BCSideID
USE MOD_Mesh_Vars     ,ONLY: nSides,nBCs,BoundaryType
USE MOD_Mesh_Vars     ,ONLY: Face_xGP
USE MOD_FV_Vars
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: t
REAL,INTENT(IN)     :: U_inner(PP_nVar,0:PP_N,0:PP_N,1:nSides)
REAL,INTENT(OUT)    :: gradU(PP_nVarPrim,0:PP_N,0:PP_N,1:nSides)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: p,q
REAL    :: U_Face_loc(PP_nVar,0:PP_N,0:PP_N)
REAL    :: UPrim_inner   (PP_nVarPrim,0:PP_N,0:PP_N)
REAL    :: UPrim_Face_loc(PP_nVarPrim,0:PP_N,0:PP_N)
INTEGER        :: iBC,iSide,SideID
INTEGER        :: BCType,BCState,nBCLoc
!==================================================================================================================================
DO iBC=1,nBCs
  IF(nBCByType(iBC).LE.0) CYCLE
  BCType =BoundaryType(iBC,BC_TYPE)
  BCState=BoundaryType(iBC,BC_STATE)
  nBCLoc =nBCByType(iBC)
 
  SELECT CASE(BCType)
  CASE(1) !Periodic already filled!
  CASE(2) ! exact BC = Dirichlet BC !!
    ! Determine the exact BC state
    DO iSide=1,nBCLoc
      SideID=BCSideID(iBC,iSide)
      DO q=0,PP_N; DO p=0,PP_N
          CALL ExactFunc(IniExactFunc,t,Face_xGP(:,p,q,1,SideID),U_Face_loc(:,p,q))
      END DO ; END DO
      gradU(:,:,:,SideID) = U_inner(:,:,:,SideID) - U_Face_loc
    END DO
    
  CASE DEFAULT ! unknown BCType
    CALL abort(__STAMP__,&
         'no BC defined in linearscalaradvection/getboundaryfvgradient.f90!',999,999.)
  END SELECT ! BCType
END DO
END SUBROUTINE GetBoundaryFVgradient
#endif


!==================================================================================================================================
!> Initialize boundary conditions
!==================================================================================================================================
SUBROUTINE FinalizeBC()
! MODULES
USE MOD_Equation_Vars,ONLY: BCData,nBCByType,BCSideID
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(BCData)
SDEALLOCATE(nBCByType)
SDEALLOCATE(BCSideID)
END SUBROUTINE FinalizeBC

END MODULE MOD_GetBoundaryFlux
