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
!TODO correct 2D implementation!
!==================================================================================================================================
!> \brief Contains the definitions of the physical fluxes of the equation system.
!>
!> The routine EvalFlux3D will compute the advection (Euler) part only, and can be called either for a single point or for
!> a volume cell. The fluxes are computed in three spatial dimension - for 2D computations, the fluxes in the third dimension
!> will always be set to 0.
!> EvalDiffFlux3D will do the same thing, but compute only the diffusive part of the fluxes. Additionally, a routine to compute
!> the fluxes on a single side is provided (used in the riemann routines).
!> The EvalEulerFlux1D routines are used in the Riemann solver, where only a flux in one spatial dimension is needed.
!>
!> The flux definitions are only done once in the single point routines, all other (side, volume) routines will simply wrap
!> to this definition.
!==================================================================================================================================
MODULE MOD_Flux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

ABSTRACT INTERFACE
  SUBROUTINE EvalFlux3DInt(NLoc,ULoc,UBaseLoc,f,g,h)
    INTEGER,INTENT(IN)                                                :: NLoc     !< Polynomial degree
    REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN)     :: ULoc     !< Solution
    REAL,DIMENSION(PP_nVarBase,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN) :: UBaseLoc !< Base flow Value
    REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT)    :: f        !< x-flux 
    REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT)    :: g        !< y-flux 
    REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT)    :: h        !< z-flux 
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  PURE SUBROUTINE EvalFlux1DInt(U,UB,F)
    REAL,INTENT(IN)     :: U(PP_nVar)   !< vector of conservative variables
    REAL,INTENT(IN)     :: UB(PP_nVarBase) !< vector of conservative variables
    REAL,INTENT(OUT)    :: F(PP_nVar)   !< Cartesian flux in "x" direction
  END SUBROUTINE
END INTERFACE

PROCEDURE(EvalFlux3DInt),POINTER :: EvalFlux3D    !< pointer defining the standard volume flux
PROCEDURE(EvalFlux1DInt),POINTER :: EvalEulerFlux1D    !< pointer defining the standard surface flux

INTEGER,PARAMETER      :: PRM_FLUX_LEE          = 1
INTEGER,PARAMETER      :: PRM_FLUX_APE          = 2

INTERFACE InitFlux
  MODULE PROCEDURE InitFlux
END INTERFACE

PUBLIC::InitFlux,EvalFlux3D,EvalEulerFlux1D,DefineParametersFlux
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersFlux()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Flux")
CALL prms%CreateIntFromStringOption('EqnType',   "System of equations used: LEE or APE",&
                                                   "LEE")
CALL addStrListEntry('EqnType','lee',               PRM_FLUX_LEE)
CALL addStrListEntry('EqnType','ape',               PRM_FLUX_APE)
END SUBROUTINE DefineParametersFlux



!==================================================================================================================================!
!> Initialize flux routines and set pointers
!==================================================================================================================================!
SUBROUTINE InitFlux()
! MODULES
USE MOD_Preproc
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                                :: EqnType !< Switch type of equation system 1: LEE, 2: APE !TODO make global?
!==================================================================================================================================
EqnType = GETINTFROMSTR('EqnType')
SELECT CASE(EqnType)
CASE(PRM_FLUX_LEE)
  EvalFlux3D      => EvalFlux3D_lee
  EvalEulerFlux1D => EvalEulerFlux1D_LEE
CASE(PRM_FLUX_APE)
  EvalFlux3D      => EvalFlux3D_ape
  EvalEulerFlux1D => EvalEulerFlux1D_APE
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
      'Type of Equation system not defined!')
END SELECT
END SUBROUTINE InitFlux


!==================================================================================================================================
!> Compute linearized Euler fluxes using the conservative variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFlux3D_LEE(NLoc,ULoc,UBaseLoc,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                       :: NLoc     !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN)  :: ULoc     ! Solution
REAL,DIMENSION(PP_nVarBase,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN) :: UBaseLoc ! Base flow
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: f        ! Cartesian flux (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: g        ! Cartesian flux (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: h        ! Cartesian flux (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: v1,v2,v3,p,rho                          ! auxiliary variables
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(NLoc);  DO j=0,NLoc; DO i=0,NLoc
  ! auxiliary variables
  rho  = ULoc(DENS,i,j,k)
  v1   = ULoc(VEL1,i,j,k)      ! u
  v2   = ULoc(VEL2,i,j,k)      ! v
  v3   = ULoc(VEL3,i,j,k)      ! w
  p    = ULoc(PRES,i,j,k)

#if PP_dim==3
  ! Euler fluxes x-direction
  f(DENS,i,j,k)=rho*UBaseLoc(VEL1,i,j,k)+UBaseLoc(DENS,i,j,k)*v1         ! u_0*rho+rho_0*u
  f(VEL1,i,j,k)=UBaseLoc(VEL1,i,j,k)*v1+(p/UBaseLoc(DENS,i,j,k))         ! u_0*u+p/rho_0
  f(VEL2,i,j,k)=UBaseLoc(VEL1,i,j,k)*v2                                  ! u_0*v
  f(VEL3,i,j,k)=UBaseLoc(VEL1,i,j,k)*v3                                  ! u_0*w
  f(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v1+UBaseLoc(VEL1,i,j,k)*p     ! gamma*p_0*u+u_0*p
  ! Euler fluxes y-direction
  g(DENS,i,j,k)=rho*UBaseLoc(VEL2,i,j,k)+UBaseLoc(DENS,i,j,k)*v2         ! v_0*rho+rho_0*v
  g(VEL1,i,j,k)=UBaseLoc(VEL2,i,j,k)*v1                                  ! v_0*u
  g(VEL2,i,j,k)=UBaseLoc(VEL2,i,j,k)*v2+(p/UBaseLoc(DENS,i,j,k))         ! v_0*v+p/rho_0
  g(VEL3,i,j,k)=UBaseLoc(VEL2,i,j,k)*v3                                  ! v_0*w
  g(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v2+UBaseLoc(VEL2,i,j,k)*p     ! gamma*p_0*v+v_0*p
  ! Euler fluxes z-direction
  h(DENS,i,j,k)=rho*UBaseLoc(VEL3,i,j,k)+UBaseLoc(DENS,i,j,k)*v3         ! w_0*rho+rho_0*w
  h(VEL1,i,j,k)=UBaseLoc(VEL3,i,j,k)*v1                                  ! w_0*u
  h(VEL2,i,j,k)=UBaseLoc(VEL3,i,j,k)*v2                                  ! w_0*v
  h(VEL3,i,j,k)=UBaseLoc(VEL3,i,j,k)*v3+(p/UBaseLoc(DENS,i,j,k))         ! w_0*w+p/rho_0
  h(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v3+UBaseLoc(VEL3,i,j,k)*p     ! gamma*p_0*w+w_0*p
#else
  ! Euler fluxes x-direction
  f(DENS,i,j,k)=rho*UBaseLoc(VEL1,i,j,k)+UBaseLoc(DENS,i,j,k)*v1         ! u_0*rho+rho_0*u
  f(VEL1,i,j,k)=UBaseLoc(VEL1,i,j,k)*v1+(p/UBaseLoc(DENS,i,j,k))         ! u_0*u+p/rho_0
  f(VEL2,i,j,k)=UBaseLoc(VEL1,i,j,k)*v2                                  ! u_0*v
  f(VEL3,i,j,k)=0.
  f(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v1+UBaseLoc(VEL1,i,j,k)*p     ! gamma*p_0*u+u_0*p
  ! Euler fluxes y-direction
  g(DENS,i,j,k)=rho*UBaseLoc(VEL2,i,j,k)+UBaseLoc(DENS,i,j,k)*v2         ! v_0*rho+rho_0*v
  g(VEL1,i,j,k)=UBaseLoc(VEL2,i,j,k)*v1                                  ! v_0*u
  g(VEL2,i,j,k)=UBaseLoc(VEL2,i,j,k)*v2+(p/UBaseLoc(DENS,i,j,k))         ! v_0*v+p/rho_0
  g(VEL3,i,j,k)=0.                                                       ! v_0*w
  g(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v2+UBaseLoc(VEL2,i,j,k)*p     ! gamma*p_0*v+v_0*p
  ! Euler fluxes z-direction
  h(:,i,j,k)=0.
#endif  
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFlux3D_LEE

!==================================================================================================================================
!> Computes 1D Euler flux
!==================================================================================================================================
PURE SUBROUTINE EvalEulerFlux1D_LEE(U,UB,F)
! MODULES
USE MOD_Equation_Vars ,ONLY:Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_nVar)   !< vector of conservative variables
REAL,INTENT(IN)     :: UB(PP_nVarBase) !< vector of conservative variables
REAL,INTENT(OUT)    :: F(PP_nVar)   !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
! Euler fluxes x-direction
F(1)= U(DENS)*UB(VEL1)+UB(DENS)*U(VEL1)         ! rho*u_0+rho_0*u
F(2)= UB(VEL1)*U(VEL1)+U(PRES)/UB(DENS)         ! u_0*u+p/rho_0
F(3)= UB(VEL1)*U(VEL2)                          ! u_0*v
#if PP_dim==3
F(4)= UB(VEL1)*U(VEL3)                          ! u_0*w
#else
F(4)=0.
#endif
F(5)= Kappa*UB(PRES)*U(VEL1)+UB(VEL1)*U(PRES)   ! gamma*p_0*u+u_0*p
END SUBROUTINE EvalEulerFlux1D_LEE


!==================================================================================================================================
!> Compute APE fluxes using the conservative variables for every volume Gauss point.
!==================================================================================================================================
SUBROUTINE EvalFlux3D_APE(NLoc,ULoc,UBaseLoc,f,g,h)
! MODULES
USE MOD_PreProc
USE MOD_Equation_Vars,ONLY:Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                       :: NLoc     !< Polynomial degree
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN)  :: ULoc     ! Solution
REAL,DIMENSION(PP_nVarBase,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(IN) :: UBaseLoc ! Base flow
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: f        ! Cartesian flux (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: g        ! Cartesian flux (iVar,i,j,k)
REAL,DIMENSION(PP_nVar,0:NLoc,0:NLoc,0:ZDIM(NLoc)),INTENT(OUT) :: h        ! Cartesian flux (iVar,i,j,k)
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: v(3),p,rho                          ! auxiliary variables
INTEGER             :: i,j,k
!==================================================================================================================================
DO k=0,ZDIM(NLoc);  DO j=0,NLoc; DO i=0,NLoc
  ! auxiliary variables
  rho  = ULoc(DENS,i,j,k)
  v(1) = ULoc(VEL1,i,j,k)      ! u
  v(2) = ULoc(VEL2,i,j,k)      ! v
  v(3) = ULoc(VEL3,i,j,k)      ! w
  p    = ULoc(PRES,i,j,k)

  ! Euler fluxes x-direction
  f(DENS,i,j,k)=0.
  f(VEL1,i,j,k)=DOT_PRODUCT(UBaseLoc(VELV,i,j,k),v)+(p/UBaseLoc(DENS,i,j,k)) ! u_0*u+v_0*v+w_0*w+p/rho_0
  f(VEL2,i,j,k)=0.
  f(VEL3,i,j,k)=0.
  f(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v(1)+UBaseLoc(VEL1,i,j,k)*p       ! gamma*p_0*u+u_0*p
  ! Euler fluxes y-direction
  g(DENS,i,j,k)=0.
  g(VEL1,i,j,k)=0.
  g(VEL2,i,j,k)=DOT_PRODUCT(UBaseLoc(VELV,i,j,k),v)+(p/UBaseLoc(DENS,i,j,k)) ! u_0*u+v_0*v+w_0*w+p/rho_0
  g(VEL3,i,j,k)=0.
  g(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v(2)+UBaseLoc(VEL2,i,j,k)*p       ! gamma*p_0*v+v_0*p
#if PP_dim==3
  ! Euler fluxes z-direction
  h(DENS,i,j,k)=0.
  h(VEL1,i,j,k)=0.
  h(VEL2,i,j,k)=0.
  h(VEL3,i,j,k)=DOT_PRODUCT(UBaseLoc(VELV,i,j,k),v)+(p/UBaseLoc(DENS,i,j,k)) ! u_0*u+v_0*v+w_0*w+p/rho_0
  h(PRES,i,j,k)=Kappa*UBaseLoc(PRES,i,j,k)*v(3)+UBaseLoc(VEL3,i,j,k)*p       ! gamma*p_0*w+w_0*p
#else
  h(:,i,j,k)=0.
#endif  
END DO; END DO; END DO ! i,j,k
END SUBROUTINE EvalFlux3D_APE

!==================================================================================================================================
!> Computes 1D APE flux
!==================================================================================================================================
PURE SUBROUTINE EvalEulerFlux1D_APE(U,UB,F)
! MODULES
USE MOD_Equation_Vars ,ONLY:Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)     :: U(PP_nVar)   !< vector of conservative variables
REAL,INTENT(IN)     :: UB(PP_nVarBase) !< vector of conservative variables
REAL,INTENT(OUT)    :: F(PP_nVar)   !< Cartesian flux in "x" direction
!----------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================================================
! Euler fluxes x-direction
F(1)= 0.
F(2)= DOT_PRODUCT(UB(VELV),U(VELV))+U(PRES)/UB(DENS)  ! u_0*u+p/rho_0
F(3)= 0.                                              ! u_0*v
F(4)= 0.
F(5)= Kappa*UB(PRES)*U(VEL1)+UB(VEL1)*U(PRES)         ! gamma*p_0*u+u_0*p
END SUBROUTINE EvalEulerFlux1D_APE

END MODULE MOD_Flux
