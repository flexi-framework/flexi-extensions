!=================================================================================================================================
! Copyright (c) 2010-2017 Prof. Claus-Dieter Munz
! Copyright (c) 2016-2017 Gregor Gassner (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Florian Hindenlang (github.com/project-fluxo/fluxo)
! Copyright (c) 2016-2017 Andrew Winters (github.com/project-fluxo/fluxo)
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

!==================================================================================================================================
!> Contains the routines for computing the fluxes for the Split-DG algorithm
!> Attention: The volume split fluxes have to incorporate a factor of 2 coming from the splitDG formulation. This can be used to
!> cancel out a single factor of 1/2 in the averages, thus reducing the required operations.
!> See the given references for details of the split fluxes or e.g. Gassner, Gregor J., Andrew R. Winters, and David A. Kopriva.
!> "Split form nodal discontinuous Galerkin schemes with summation-by-parts property for the compressible Euler equations."
!> Journal of Computational Physics 327 (2016): 39-66. for an overview.
!==================================================================================================================================
#include "flexi.h"
#include "eos.h"

MODULE MOD_SplitFlux
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PPURE SUBROUTINE VolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
    REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef,U
    REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef,UPrim
    REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef,MeshVel
    REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef,M
    REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  PPURE SUBROUTINE SurfaceFlux(U_LL,U_RR,MeshVel,F)
    REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL,U_RR
    REAL,                   INTENT(IN)  :: MeshVel
    REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(VolumeFlux),POINTER    :: SplitDGVolume_pointer    !< pointer defining the SpliDG formulation beeing used
PROCEDURE(SurfaceFlux),POINTER   :: SplitDGSurface_pointer   !< pointer defining the SpliDG formulation beeing used

INTEGER,PARAMETER      :: PRM_SPLITDG_SD          = 0
INTEGER,PARAMETER      :: PRM_SPLITDG_MO          = 1
INTEGER,PARAMETER      :: PRM_SPLITDG_DU          = 2
INTEGER,PARAMETER      :: PRM_SPLITDG_KG          = 3
INTEGER,PARAMETER      :: PRM_SPLITDG_PI          = 4
INTEGER,PARAMETER      :: PRM_SPLITDG_CH          = 5

INTERFACE InitSplitDG
  MODULE PROCEDURE InitSplitDG
END INTERFACE

PUBLIC::InitSplitDG,DefineParametersSplitDG
PUBLIC::SplitDGSurface_pointer,SplitDGVolume_pointer
PUBLIC::GetLogMean
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersSplitDG()
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
CALL prms%SetSection("SplitDG")
CALL prms%CreateIntFromStringOption('SplitDG',"SplitDG formulation to be used: SD, MO, DU, KG, PI, CH","PI")
CALL addStrListEntry('SplitDG','sd',           PRM_SPLITDG_SD)
CALL addStrListEntry('SplitDG','mo',           PRM_SPLITDG_MO)
CALL addStrListEntry('SplitDG','du',           PRM_SPLITDG_DU)
CALL addStrListEntry('SplitDG','kg',           PRM_SPLITDG_KG)
CALL addStrListEntry('SplitDG','pi',           PRM_SPLITDG_PI)
CALL addStrListEntry('SplitDG','ch',           PRM_SPLITDG_CH)

END SUBROUTINE DefineParametersSplitDG

!==================================================================================================================================!
!> Initialize function pointers for the specific split version in use
!==================================================================================================================================!
SUBROUTINE InitSplitDG()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                     :: SplitDG
!==================================================================================================================================
! check if Gauss-Lobatto-Pointset is beeing used
#if (PP_NodeType==1)
CALL CollectiveStop(__STAMP__,&
  'Wrong Pointset: Gauss-Lobatto-Points are mandatory for using SplitDG !')
#endif
! set pointers
SplitDG = GETINTFROMSTR('SplitDG')
SELECT CASE(SplitDG)
CASE(PRM_SPLITDG_SD)
  SplitDGVolume_pointer  => SplitVolumeFluxSD
  SplitDGSurface_pointer => SplitSurfaceFluxSD
CASE(PRM_SPLITDG_MO)
  SplitDGVolume_pointer  => SplitVolumeFluxMO
  SplitDGSurface_pointer => SplitSurfaceFluxMO
CASE(PRM_SPLITDG_DU)
  SplitDGVolume_pointer  => SplitVolumeFluxDU
  SplitDGSurface_pointer => SplitSurfaceFluxDU
CASE(PRM_SPLITDG_KG)
  SplitDGVolume_pointer  => SplitVolumeFluxKG
  SplitDGSurface_pointer => SplitSurfaceFluxKG
CASE(PRM_SPLITDG_PI)
  SplitDGVolume_pointer  => SplitVolumeFluxPI
  SplitDGSurface_pointer => SplitSurfaceFluxPI
CASE(PRM_SPLITDG_CH)
  SplitDGVolume_pointer  => SplitVolumeFluxCH
  SplitDGSurface_pointer => SplitSurfaceFluxCH
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'SplitDG formulation not defined!')
END SELECT
END SUBROUTINE InitSplitDG

!==================================================================================================================================
!> Computes the Split-Flux retaining the standard NS-Equations
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxSD(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel       !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: rhoEpRef,rhoEp ! auxiliary variable for (rho*E+p)
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde  ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)                 :: hTilde         ! flux in physical space
#endif
!==================================================================================================================================
! compute auxiliary variables, total energy density plus pressure
rhoEpRef = URef(5) + UPrimRef(5)
rhoEp    = U(5) + UPrim(5)

! local Euler fluxes x-direction
fTilde(1) = (URef(2) + U(2))                                           ! {rho*u}
fTilde(2) = (URef(2)*UPrimRef(2)+UPrimRef(5) + U(2)*UPrim(2)+UPrim(5)) ! {rho*u²}+{p}
fTilde(3) = (URef(2)*UPrimRef(3) + U(2)*UPrim(3))                      ! {rho*u*v}
#if PP_dim == 3
fTilde(4) = (URef(2)*UPrimRef(4) + U(2)*UPrim(4))                      ! {rho*u*w}
#else
fTilde(4) = 0.
#endif
fTilde(5) = (rhoEpRef*UPrimRef(2) + rhoEp*UPrim(2))                    ! {(rho*E+p)*u}
! local Euler fluxes y-direction
gTilde(1) = (URef(3) + U(3))                                           ! {rho*v}
gTilde(2) = (URef(2)*UPrimRef(3) + U(2)*UPrim(3))                      ! {rho*u*v}
gTilde(3) = (URef(3)*UPrimRef(3)+UPrimRef(5) + U(3)*UPrim(3)+UPrim(5)) ! {rho*v²}+{p}
#if PP_dim == 3
gTilde(4) = (URef(3)*UPrimRef(4) + U(3)*UPrim(4))                      ! {rho*v*w}
#else
gTilde(4) = 0.
#endif
gTilde(5) = (rhoEpRef*UPrimRef(3) + rhoEp*UPrim(3))                    ! {(rho*E+p)*v}
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) = (URef(4) + U(4))                                           ! {rho*w}
hTilde(2) = (URef(2)*UPrimRef(4) + U(2)*UPrim(4))                      ! {rho*u*w}
hTilde(3) = (URef(3)*UPrimRef(4) + U(3)*UPrim(4))                      ! {rho*v*w}
hTilde(4) = (URef(4)*UPrimRef(4)+UPrimRef(5) + U(4)*UPrim(4)+UPrim(5)) ! {rho*v²+p}
hTilde(5) = (rhoEpRef*UPrimRef(4) + rhoEp*UPrim(4))                    ! {(rho*E+p)*w}
#endif

! ALE fluxes
fTilde(1) = fTilde(1)-0.5*(MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))
fTilde(2) = fTilde(2)-0.5*(MeshVelRef(1)+MeshVel(1))*(URef(2)+U(2))
fTilde(3) = fTilde(3)-0.5*(MeshVelRef(1)+MeshVel(1))*(URef(3)+U(3))
#if PP_dim == 3
fTilde(4) = fTilde(4)-0.5*(MeshVelRef(1)+MeshVel(1))*(URef(4)+U(4))
#endif
fTilde(5) = fTilde(5)-0.5*(MeshVelRef(1)+MeshVel(1))*(URef(5)+U(5))

gTilde(1) = gTilde(1)-0.5*(MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))
gTilde(2) = gTilde(2)-0.5*(MeshVelRef(2)+MeshVel(2))*(URef(2)+U(2))
gTilde(3) = gTilde(3)-0.5*(MeshVelRef(2)+MeshVel(2))*(URef(3)+U(3))
#if PP_dim == 3
gTilde(4) = gTilde(4)-0.5*(MeshVelRef(2)+MeshVel(2))*(URef(4)+U(4))
#endif
gTilde(5) = gTilde(5)-0.5*(MeshVelRef(2)+MeshVel(2))*(URef(5)+U(5))

#if PP_dim == 3
hTilde(1) = hTilde(1)-0.5*(MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))
hTilde(2) = hTilde(2)-0.5*(MeshVelRef(3)+MeshVel(3))*(URef(2)+U(2))
hTilde(3) = hTilde(3)-0.5*(MeshVelRef(3)+MeshVel(3))*(URef(3)+U(3))
hTilde(4) = hTilde(4)-0.5*(MeshVelRef(3)+MeshVel(3))*(URef(4)+U(4))
hTilde(5) = hTilde(5)-0.5*(MeshVelRef(3)+MeshVel(3))*(URef(5)+U(5))
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxSD

!==================================================================================================================================
!> Computes the surface flux for the split formulation retaining the standard NS-Equations
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxSD(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================

F(1)= 0.5*(U_LL(MOM1)+U_RR(MOM1))                                                ! {rho*u}
F(2)= 0.5*(U_LL(MOM1)*U_LL(VEL1)+U_LL(PRES)+U_RR(MOM1)*U_RR(VEL1)+U_RR(PRES))    ! {rho*u²}+{p}
F(3)= 0.5*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))                          ! {rho*u*v}
#if PP_dim == 3
F(4)= 0.5*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))                          ! {rho*u*w}
#else
F(4)= 0.
#endif
F(5)= 0.5*((U_LL(ENER)+U_LL(PRES))*U_LL(VEL1)+(U_RR(ENER)+U_RR(PRES))*U_RR(VEL1))! {(rho*E+p)*u}


! ALE fluxes
F(1)=F(1)-0.5*MeshVel*(U_LL(DENS)+U_RR(DENS))
F(2)=F(2)-0.5*MeshVel*(U_LL(MOM1)+U_RR(MOM1))
F(3)=F(3)-0.5*MeshVel*(U_LL(MOM2)+U_RR(MOM2))
#if PP_dim == 3
F(4)=F(4)-0.5*MeshVel*(U_LL(MOM3)+U_RR(MOM3))
#endif
F(5)=F(5)-0.5*MeshVel*(U_LL(ENER)+U_RR(ENER))
END SUBROUTINE SplitSurfaceFluxSD

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Ducros
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Uses quadratic split forms in all equations. Reference: Ducros, F., et al. "High-order fluxes for conservative
!> skew-symmetric-like schemes in structured meshes: application to compressible flows." 
!> Journal of Computational Physics 161.1 (2000): 114-139.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxDU(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel       !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(PP_nVar)             :: fTilde,gTilde     ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)             :: hTilde            ! flux in physical space
#endif
!==================================================================================================================================

! local Euler fluxes x-direction
fTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                          ! {rho}*{u}
fTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(2)+UPrim(2)) + (UPrimRef(5)+UPrim(5)) ! {rho*u}*{u}+{p}
fTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(2)+UPrim(2))                          ! {rho*v}*{u}
#if PP_dim == 3
fTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(2)+UPrim(2))                          ! {rho*w}*{u}
#else
fTilde(4) = 0.
#endif
fTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(2)+UPrim(2))     ! ({rho*E}+{p})*{u}
! local Euler fluxes y-direction
gTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                          ! {rho}*{v}
gTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(3)+UPrim(3))                          ! {rho*u}*{v}
gTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(3)+UPrim(3)) + (UPrimRef(5)+UPrim(5)) ! {rho*v}*{v}+{p}
#if PP_dim == 3
gTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(3)+UPrim(3))                          ! {rho*w}*{v}
#else
gTilde(4) = 0.
#endif
gTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(3)+UPrim(3))     ! ({rho*E}+{p})*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) = 0.5*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                          ! {rho}*{w}
hTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(4)+UPrim(4))                          ! {rho*u}*{w}
hTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(4)+UPrim(4))                          ! {rho*v}*{w}
hTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(4)+UPrim(4)) + (UPrimRef(5)+UPrim(5)) ! {rho*w}*{w}+{p}
hTilde(5) = 0.5*(URef(5)+U(5)+UPrimRef(5)+UPrim(5))*(UPrimRef(4)+UPrim(4))     ! ({rho*E}+{p})*{w}
#endif

! Kinetic energy stable ALE fluxes
CALL KepALEVolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,fTilde,gTilde&
#if PP_dim == 3
                                  ,hTilde&
#endif
                                 )

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxDU

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Ducros
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxDU(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
F(1)= 0.25*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho}*{u}
F(2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(3)= 0.25*(U_LL(MOM2)+U_RR(MOM2))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*v}*{u}
#if PP_dim == 3
F(4)= 0.25*(U_LL(MOM3)+U_RR(MOM3))*(U_LL(VEL1)+U_RR(VEL1))                               ! {rho*w}*{u}
#else
F(4)= 0.
#endif
F(5)= 0.25*(U_LL(ENER)+U_RR(ENER)+U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))         ! ({rho*E}+{p})*{u}

! Kinetic energy stable ALE fluxes
CALL KepALESurfaceFlux(U_LL,U_RR,MeshVel,F)

END SUBROUTINE SplitSurfaceFluxDU

!==================================================================================================================================
!> Computes the Split-Flux retaining the KEP formulation of Kennedy and Gruber
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Reference: Kennedy, Christopher A., and Andrea Gruber. "Reduced aliasing formulations of the convective terms within the 
!> Navier–Stokes equations for a compressible fluid." Journal of Computational Physics 227.3 (2008): 1676-1700.
!> Uses a quadratic splitting for u*p and a cubic splitting for rho*e*u in the energy equation. 
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxKG(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel       !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: E,ERef        ! auxiliary variables for the specific total energy
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)                 :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! specific total energy
ERef = URef(5)/URef(1)
E    = U(5)/U(1)

! local Euler fluxes x-direction
fTilde(1) = 0.5* (URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                             ! {rho}*{u}
fTilde(2) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{u}²+{p}
fTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(3)+UPrim(3))      ! {rho}*{u}*{v}
#if PP_dim == 3
fTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(4)+UPrim(4))      ! {rho}*{u}*{w}
#else
fTilde(4) = 0.
#endif
fTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(eRef+e) + &
            0.5* (UPrimRef(5)+UPrim(5))*(UPrimRef(2)+UPrim(2))                     ! {rho}*{E}*{u}+{p}*{u}
! local Euler fluxes y-direction
gTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                             ! {rho}*{v}
gTilde(2) = fTilde(3)                                                              ! {rho}*{v}*{u}
gTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{v}²+{p}
#if PP_dim == 3
gTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(UPrimRef(4)+UPrim(4))      ! {rho}*{v}*{w}
#else
gTilde(4) = 0.
#endif
gTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(eRef+e) + &
            0.5* (UPrimRef(5)+UPrim(5))*(UPrimRef(3)+UPrim(3))                     ! {rho}*{E}*{v}+{p}*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                             ! {rho}*{w}
hTilde(2) = fTilde(4)                                                              ! {rho}*{w}*{u}
hTilde(3) = gTilde(4)                                                              ! {rho}*{w}*{v}
hTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{w}²+{p}
hTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))*(eRef+e) + &
            0.5 *(UPrimRef(5)+UPrim(5))*(UPrimRef(4)+UPrim(4))                     ! {rho}*{E}*{w}+{p}*{w}
#endif

! Kinetic energy stable ALE fluxes
CALL KepALEVolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,fTilde,gTilde&
#if PP_dim == 3
                                  ,hTilde&
#endif
                                 )
! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxKG

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Kennedy and Gruber
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxKG(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: E_LL,E_RR ! auxiliary variables for the specific total energy
!==================================================================================================================================
! specific total energy
E_LL = U_LL(ENER)/U_LL(DENS)
E_RR = U_RR(ENER)/U_RR(DENS)
!compute flux
F(1)= 0.25* (U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{u}
F(2)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))**2 + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(3)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL2)+U_RR(VEL2))          ! {rho}*{u}*{v}
#if PP_dim == 3
F(4)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL3)+U_RR(VEL3))          ! {rho}*{u}*{w}
#else
F(4)= 0.
#endif
F(5)= 0.125*(U_LL(DENS)+U_RR(DENS))*(E_LL+E_RR)*(U_LL(VEL1)+U_RR(VEL1)) + &
      0.25 *(U_LL(PRES)+U_RR(PRES))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{E}*{u}+{p}*{u}

! Kinetic energy stable ALE fluxes
CALL KepALESurfaceFlux(U_LL,U_RR,MeshVel,F)

END SUBROUTINE SplitSurfaceFluxKG

!==================================================================================================================================
!> Computes the Split-Flux retaining the formulation of Morinishi
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> Uses quadratic splittings in the momentum equations, but pure advection form in the energy equation. This seems to lead to a
!> rather unstable formulation!
!> Reference: Morinishi, Yohei. "Skew-symmetric form of convective terms and fully conservative finite difference schemes for
!> variable density low-Mach number flows." Journal of Computational Physics 229.2 (2010): 276-300.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxMO(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel       !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: rhoepRef,rhoep ! auxiliary variable for rho*inner energy + pressure
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde  ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)                 :: hTilde         ! flux in physical space
#endif
!==================================================================================================================================

! rho*internal energy + pressure
rhoepRef = URef(5)-0.5*URef(1)*(UPrimRef(2)**2+UPrimRef(3)**2+UPrimRef(4)**2)+UPrimRef(5)
rhoep    = U(5)-0.5*U(1)*(UPrim(2)**2+UPrim(3)**2+UPrim(4)**2)+UPrim(5)

! local Euler fluxes x-direction
fTilde(1) =     (URef(2)+U(2))                                                 ! {rho*u}
fTilde(2) = 0.5*(URef(2)+U(2))*(UPrimRef(2)+UPrim(2)) + (UPrimRef(5)+UPrim(5)) ! {rho*u}*{u}+{p}
fTilde(3) = 0.5*(Uref(2)+U(2))*(UPrimRef(3)+UPrim(3))                          ! {rho*u}*{v}
#if PP_dim == 3
fTilde(4) = 0.5*(URef(2)+U(2))*(UPrimRef(4)+UPrim(4))                          ! {rho*u}*{w}
#else
fTilde(4) = 0.
#endif
fTilde(5) = (rhoepRef*UPrimRef(2)+rhoep*UPrim(2)) + &                          !{(rho*e+p)*u} +
            0.5*(URef(2)*UPrimRef(2)+U(2)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*u²}*{u} +
            0.5*(URef(2)*UPrimRef(3)+U(2)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*u*v}*{v} +
            0.5*(URef(2)*UPrimRef(4)+U(2)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*u*w}*{w} -
            0.5*(URef(2)*UPrimRef(2)*UPrimRef(2)+U(2)*UPrim(2)*UPrim(2))   - & !1/2*({rho*u³} +
            0.5*(URef(2)*UPrimRef(3)*UPrimRef(3)+U(2)*UPrim(3)*UPrim(3))   - & !{rho*u*v²} +
            0.5*(URef(2)*UPrimRef(4)*UPrimRef(4)+U(2)*UPrim(4)*UPrim(4))       !{rho*u*w²})
! local Euler fluxes y-direction
gTilde(1) =     (URef(3)+U(3))                                                 ! {rho*v}
gTilde(2) = 0.5*(Uref(3)+U(3))*(UPrimRef(2)+UPrim(2))                          ! {rho*v}*{u}
gTilde(3) = 0.5*(URef(3)+U(3))*(UPrimRef(3)+UPrim(3)) + (UPrimRef(5)+UPrim(5)) ! {rho*v}*{v}+{p}
#if PP_dim == 3
gTilde(4) = 0.5*(URef(3)+U(3))*(UPrimRef(4)+UPrim(4))                          ! {rho*v}*{w}
#else
gTilde(4) = 0.
#endif
gTilde(5) = (rhoepRef*UPrimRef(3)+rhoep*UPrim(3)) + &                          !{(rho*e+p)*v} +
            0.5*(URef(3)*UPrimRef(2)+U(3)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*v*u}*{u} +
            0.5*(URef(3)*UPrimRef(3)+U(3)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*v²}*{v} +
            0.5*(URef(3)*UPrimRef(4)+U(3)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*v*w}*{w} -
            0.5*(URef(3)*UPrimRef(2)*UPrimRef(2)+U(3)*UPrim(2)*UPrim(2))   - & !1/2*({rho*v*u²} +
            0.5*(URef(3)*UPrimRef(3)*UPrimRef(3)+U(3)*UPrim(3)*UPrim(3))   - & !{rho*v³} +
            0.5*(URef(3)*UPrimRef(4)*UPrimRef(4)+U(3)*UPrim(4)*UPrim(4))       !{rho*v*w²})
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) =     (URef(4)+U(4))                                                 ! {rho*w}
hTilde(2) = 0.5*(Uref(4)+U(4))*(UPrimRef(2)+UPrim(2))                          ! {rho*w}*{u}
hTilde(3) = 0.5*(URef(4)+U(4))*(UPrimRef(3)+UPrim(3))                          ! {rho*w}*{v}
hTilde(4) = 0.5*(URef(4)+U(4))*(UPrimRef(4)+UPrim(4)) + (UPrimRef(5)+UPrim(5)) ! {rho*w}*{w}+{p}
hTilde(5) = (rhoepRef*UPrimRef(4)+rhoep*UPrim(4)) + &                          !{(rho*e+p)*w} +
            0.5*(URef(4)*UPrimRef(2)+U(4)*UPrim(2))*(UPrimRef(2)+UPrim(2)) + & !{rho*w*u}*{u} +
            0.5*(URef(4)*UPrimRef(3)+U(4)*UPrim(3))*(UPrimRef(3)+UPrim(3)) + & !{rho*w*v}*{v} +
            0.5*(URef(4)*UPrimRef(4)+U(4)*UPrim(4))*(UPrimRef(4)+UPrim(4)) - & !{rho*w²}*{w} -
            0.5*(URef(4)*UPrimRef(2)*UPrimRef(2)+U(4)*UPrim(2)*UPrim(2))   - & !1/2*({rho*w*u²} +
            0.5*(URef(4)*UPrimRef(3)*UPrimRef(3)+U(4)*UPrim(3)*UPrim(3))   - & !{rho*w*v²} +
            0.5*(URef(4)*UPrimRef(4)*UPrimRef(4)+U(4)*UPrim(4)*UPrim(4))       !{rho*w³})
#endif

! Kinetic energy stable ALE fluxes
CALL KepALEVolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,fTilde,gTilde&
#if PP_dim == 3
                                  ,hTilde&
#endif
                                 )

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxMO

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Morinishi
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxMO(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F    !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: rhoep_LL,rhoep_RR
!==================================================================================================================================
! rho*internal energy + pressure
rhoep_LL = U_LL(ENER)-0.5*U_LL(DENS)*(U_LL(VEL1)**2+U_LL(VEL2)**2+U_LL(VEL3)**2)+U_LL(PRES)
rhoep_RR = U_RR(ENER)-0.5*U_RR(DENS)*(U_RR(VEL1)**2+U_RR(VEL2)**2+U_RR(VEL3)**2)+U_RR(PRES)

! compute flux
F(1)= 0.5 *(U_LL(MOM1)+U_RR(MOM1))                                                       ! {rho*u}
F(2)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL1)+U_RR(VEL1)) + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho*u}*{u}+{p}
F(3)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL2)+U_RR(VEL2))                               ! {rho*u}*{v}
#if PP_dim == 3
F(4)= 0.25*(U_LL(MOM1)+U_RR(MOM1))*(U_LL(VEL3)+U_RR(VEL3))                               ! {rho*u}*{w}
#else
F(4)= 0.
#endif
F(5)= 0.5 *(rhoep_LL*U_LL(VEL1)+rhoep_RR*U_RR(VEL1)) +  &                                !{(rho*e+p)*u} +
      0.25*(U_LL(MOM1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1))*(U_LL(VEL1)+U_RR(VEL1)) + &     !{rho*u²}*{u} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2))*(U_LL(VEL2)+U_RR(VEL2)) + &     !{rho*u*v}*{v} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3))*(U_LL(VEL3)+U_RR(VEL3)) - &     !{rho*u*w}*{w} -
      0.25*(U_LL(MOM1)*U_LL(VEL1)*U_LL(VEL1)+U_RR(MOM1)*U_RR(VEL1)*U_RR(VEL1)) - &       !1/2*({rho*u³} +
      0.25*(U_LL(MOM1)*U_LL(VEL2)*U_LL(VEL2)+U_RR(MOM1)*U_RR(VEL2)*U_RR(VEL2)) - &       !{rho*u*v²} +
      0.25*(U_LL(MOM1)*U_LL(VEL3)*U_LL(VEL3)+U_RR(MOM1)*U_RR(VEL3)*U_RR(VEL3))           !{rho*u*w²})

! Kinetic energy stable ALE fluxes
CALL KepALESurfaceFlux(U_LL,U_RR,MeshVel,F)

END SUBROUTINE SplitSurfaceFluxMO

!==================================================================================================================================
!> Computes the Split-Flux retaining the KEP formulation of Pirozzoli
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> This is a slight variation of the form proposed by Kennedy and Gruber, and is based on a cubic splitting of rho*H*u in the energy
!> equation. It was presented in the (worth reading) overview: Pirozzoli, Sergio. "Numerical methods for high-speed flows."
!> Annual review of fluid mechanics 43 (2011): 163-194.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxPI(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef          !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef      !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim         !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel       !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef          !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M             !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux          !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: H,HRef        ! auxiliary variables for the specific enthalpy
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)                 :: hTilde        ! flux in physical space
#endif
!==================================================================================================================================

! specific enthalpy, H=E+p/rho=(rhoE+p)/rho
HRef = (URef(5)+UPrimRef(5))/URef(1)
H    = (U(5)+UPrim(5))/U(1)

! local Euler fluxes x-direction
fTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))                             ! {rho}*{u}
fTilde(2) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{u}²+{p}
fTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(3)+UPrim(3))      ! {rho}*{u}*{v}
#if PP_dim == 3
fTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(UPrimRef(4)+UPrim(4))      ! {rho}*{u}*{w}
#else
fTilde(4) = 0.
#endif
fTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))*(HRef+H)                    ! {rho}*{H}*{u}
! local Euler fluxes y-direction
gTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))                             ! {rho}*{v}
gTilde(2) = fTilde(3)                                                              ! {rho}*{v}*{u}
gTilde(3) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{v}²+{p}
#if PP_dim == 3
gTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(UPrimRef(4)+UPrim(4))      ! {rho}*{v}*{w}
#else
gTilde(4) = 0.
#endif
gTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))*(HRef+H)                    ! {rho}*{H}*{v}
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) = 0.5 *(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))                             ! {rho}*{w}
hTilde(2) = fTilde(4)                                                              ! {rho}*{w}*{u}
hTilde(3) = gTilde(4)                                                              ! {rho}*{w}*{v}
hTilde(4) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))**2 + (UPrimRef(5)+UPrim(5)) ! {rho}*{w}²+{p}
hTilde(5) = 0.25*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))*(HRef+H)                    ! {rho}*{H}*{w}
#endif

! Kinetic energy stable ALE fluxes
CALL KepALEVolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,fTilde,gTilde&
#if PP_dim == 3
                                  ,hTilde&
#endif
                                 )

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

END SUBROUTINE SplitVolumeFluxPI

!==================================================================================================================================
!> Computes the surface flux for the split formulation of Pirozzoli
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxPI(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL      !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR      !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F         !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: H_LL,H_RR ! auxiliary variables for the specific energy
!==================================================================================================================================
! specific energy, H=E+p/rho=(rhoE+p)/rho
H_LL = (U_LL(ENER)+U_LL(PRES))/U_LL(DENS)
H_RR = (U_RR(ENER)+U_RR(PRES))/U_RR(DENS)
!compute flux
F(1)= 0.25* (U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))                                  ! {rho}*{u}
F(2)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))**2 + 0.5*(U_LL(PRES)+U_RR(PRES)) ! {rho}*{u}²+{p}
F(3)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL2)+U_RR(VEL2))          ! {rho}*{u}*{v}
#if PP_dim == 3
F(4)= 0.125*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))*(U_LL(VEL3)+U_RR(VEL3))          ! {rho}*{u}*{w}
#else
F(4)= 0.
#endif
F(5)= 0.125*(U_LL(DENS)+U_RR(DENS))*(H_LL+H_RR)*(U_LL(VEL1)+U_RR(VEL1))                      ! {rho}*{H}*{u}

! Kinetic energy stable ALE fluxes
CALL KepALESurfaceFlux(U_LL,U_RR,MeshVel,F)

END SUBROUTINE SplitSurfaceFluxPI

!==================================================================================================================================
!> Computes the Split-Flux retaining the entropy conserving (and formally KEP) formulation of Chandrashekar
!> Attention 1: Factor 2 from differentiation matrix is already been considered
!> The flux after Chanrashekar uses a special computation of the pressure, based on the averages of density and inverse 
!> temperature, which correspondonds to using the harmonic average of the temperature when applying the ideal gas law.
!> Reference: Chandrashekar, Praveen. "Kinetic energy preserving and entropy stable finite volume schemes for compressible Euler
!> and Navier-Stokes equations." Communications in Computational Physics 14.5 (2013): 1252-1286.
!==================================================================================================================================
PPURE SUBROUTINE SplitVolumeFluxCH(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,MRef,M,Flux)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars, ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: URef       !< conserved variables
REAL,DIMENSION(PP_nVar    ),INTENT(IN)  :: U          !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrimRef   !< primitive variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)  :: UPrim      !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVelRef !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MeshVel    !< mesh velocities
REAL,DIMENSION(1:3        ),INTENT(IN)  :: MRef       !< metric terms
REAL,DIMENSION(1:3        ),INTENT(IN)  :: M          !< metric terms
REAL,DIMENSION(PP_nVar    ),INTENT(OUT) :: Flux       !< flux in reverence space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                    :: beta,betaRef            ! auxiliary variables for the inverse Temperature
REAL                                    :: pHatMean,HMeanMod          ! auxiliary variable for the mean pressure and specific enthalpy
REAL                                    :: uMean,vMean,wMean       ! auxiliary variable for the average velocities
REAL                                    :: uRelMean,vRelMean,wRelMean ! auxiliary variable for the average relative velocities
REAL                                    :: rhoLogMean,betaLogMean  ! auxiliary variable for the logarithmic means
REAL                                    :: rhoMean                 ! auxiliary variable for the mean density
REAL,DIMENSION(PP_nVar)                 :: fTilde,gTilde           ! flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)                 :: hTilde                  ! flux in physical space
#endif
!==================================================================================================================================
! average velocities
uMean = 0.5*(UPrimRef(2) + UPrim(2))
vMean = 0.5*(UPrimRef(3) + UPrim(3))
wMean = 0.5*(UPrimRef(4) + UPrim(4))
uRelMean = 0.5*(UPrimRef(2)-MeshVelRef(1) + UPrim(2)-MeshVel(1))
vRelMean = 0.5*(UPrimRef(3)-MeshVelRef(2) + UPrim(3)-MeshVel(2))
wRelMean = 0.5*(UPrimRef(4)-MeshVelRef(3) + UPrim(4)-MeshVel(3))

! inverse temperature
betaRef  = 0.5*URef(1)/UPrimRef(5)
beta     = 0.5*U(1)/UPrim(5)

! Density and inverse temperature logarithmic average
CALL GetLogMean(URef(1),U(1),rhoLogMean)
CALL GetLogMean(betaRef,beta,betaLogMean)
! Average of pressure and enthalpy - enthalpy is modified so we can use it with the mesh velocity
pHatMean = 0.5*(URef(1)+U(1))/(betaRef+beta)
HMeanMod = 0.5*sKappaM1/betaLogMean + &
           0.5*(UPrimRef(2)*UPrim(2) + UPrimRef(3)*UPrim(3) + UPrimRef(4)*UPrim(4))

! local Euler fluxes x-direction
fTilde(1) = rhoLogMean*uRelMean                                                    ! {rho}_log*{u-u_m}
fTilde(2) = rhoLogMean*uRelMean*uMean + pHatMean                                      ! {rho}_log*{u-u_m}*{u}+{p}
fTilde(3) = rhoLogMean*uRelMean*vMean                                              ! {rho}_log*{u-u_m}*{v}
#if PP_dim == 3
fTilde(4) = rhoLogMean*uRelMean*wMean                                              ! {rho}_log*{u_m}*{w}
#else
fTilde(4) = 0.
#endif
! {rho}_log*{u-u_m}/(2*(kappa-1*{beta}_log))+0.5*{rho}_log*{u-u_m}*{|u|^2}+{rho}*{u}/(2*{beta})
fTilde(5) = rhoLogMean*HMeanMod*uRelMean + pHatMean*uMean
! local Euler fluxes y-direction
gTilde(1) = rhoLogMean*vRelMean                                                    ! {rho}_log*{v-v_m}
gTilde(2) = rhoLogMean*vRelMean*uMean                                              ! {rho}_log*{v-v_m}*{u}
gTilde(3) = rhoLogMean*vRelMean*vMean + pHatMean                                      ! {rho}_log*{v-v_m}*{v_m}+{p}
#if PP_dim == 3
gTilde(4) = rhoLogMean*vRelMean*wMean                                              ! {rho}_log*{v-v_m}*{w}
#else
gTilde(4) = 0.
#endif
! {rho}_log*{v-v_m}/(2*(kappa-1*{beta}_log))+0.5*{rho}_log*{v-v_m}*{|u|^2}+{rho}*{v}/(2*{beta})
gTilde(5) = rhoLogMean*HMeanMod*vRelMean + pHatMean*vMean
#if PP_dim == 3
! local Euler fluxes z-direction
hTilde(1) = rhoLogMean*wRelMean                                                       ! {rho}_log*{w-w_m}
hTilde(2) = rhoLogMean*wRelMean*uMean                                                 ! {rho}_log*{w-w_m}*{u}
hTilde(3) = rhoLogMean*wRelMean*vMean                                                 ! {rho}_log*{w-w_m}*{v}
hTilde(4) = rhoLogMean*wRelMean*wMean + pHatMean                                         ! {rho}_log*{w-w_m}*{w}+{p}
! {rho}_log*{w-w_m}/(2*(kappa-1*{beta}_log))+0.5*{rho}_log*{w-w_m}*{|u|^2}+{rho}*{w}/(2*{beta})
hTilde(5) = rhoLogMean*HMeanMod*wRelMean + pHatMean*wMean
#endif

! transform into reference space
Flux(:) = 0.5*(MRef(1)+M(1))*fTilde(:) + &
#if PP_dim == 3
          0.5*(MRef(3)+M(3))*hTilde(:) + &
#endif
          0.5*(MRef(2)+M(2))*gTilde(:)

! Acount for factor of 2
Flux = Flux*2.

END SUBROUTINE SplitVolumeFluxCH

!==================================================================================================================================
!> Computes the surface flux for the entropy conserving formulation of Chandrashekar.
!> The flux after Chanrashekar uses a special computation of the pressure, based on the averages of density and inverse 
!> temperature, which correspondonds to using the harmonic average of the temperature when applying the ideal gas law.
!==================================================================================================================================
PPURE SUBROUTINE SplitSurfaceFluxCH(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars, ONLY:sKappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_LL    !< variables at the left surfaces
REAL,DIMENSION(PP_2Var),INTENT(IN)  :: U_RR    !< variables at the right surfaces
REAL,                   INTENT(IN)  :: MeshVel !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(OUT) :: F       !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                                :: beta_LL,beta_RR        ! auxiliary variables for the inverse Temperature
REAL                                :: pHatMean,HMeanMod         ! auxiliary variable for the mean pressure and enthalpy
REAL                                :: uMean,vMean,wMean,uRelMean ! auxiliary variable for the average velocities
REAL                                :: rhoLogMean,betaLogMean ! auxiliary variable for the logarithmic mean
!==================================================================================================================================
! average velocities
uMean = 0.5*(U_LL(VEL1) + U_RR(VEL1))
vMean = 0.5*(U_LL(VEL2) + U_RR(VEL2))
wMean = 0.5*(U_LL(VEL3) + U_RR(VEL3))
uRelMean = 0.5*(U_LL(VEL1) - MeshVel + U_RR(VEL1) - MeshVel)

! inverse temperature
beta_LL = 0.5*U_LL(DENS)/U_LL(PRES)
beta_RR = 0.5*U_RR(DENS)/U_RR(PRES)

! average pressure, enthalpy, density and inverse temperature
! logarithmic mean
CALL GetLogMean(U_LL(DENS),U_RR(DENS),rhoLogMean)
CALL GetLogMean(beta_LL,beta_RR,betaLogMean)
! "standard" average
pHatMean = 0.5*(U_LL(DENS)+U_RR(DENS))/(beta_LL+beta_RR)
HMeanMod = 0.5*sKappaM1/betaLogMean + &
           0.5*(U_LL(VEL1)*U_RR(VEL1) + U_LL(VEL2)*U_RR(VEL2) + U_LL(VEL3)*U_RR(VEL3))

!compute flux
F(1) = rhoLogMean*uRelMean
F(2) = F(1)*uMean + pHatMean
F(3) = F(1)*vMean
F(4) = F(1)*wMean
F(5) = F(1)*HMeanMod + pHatMean * uMean

END SUBROUTINE SplitSurfaceFluxCH

!==================================================================================================================================
!> auxilary function for calculating the logarithmic mean numerically stable according to Ismail and Roe
!==================================================================================================================================
ELEMENTAL SUBROUTINE GetLogMean(U_L,U_R,UMean)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)  :: U_L   !< variables at the left surfaces
REAL,INTENT(IN)  :: U_R   !< variables at the right surfaces
REAL,INTENT(OUT) :: UMean !< resulting flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,PARAMETER   :: epsilon = 0.01
REAL             :: chi,f,u,N ! auxiliary variables
!==================================================================================================================================
chi = U_L/U_R
f = (chi-1)/(chi+1)
u = f*f

IF (u .LT. epsilon) THEN
  N = 1.0+u/3.0+u*u/5.0+u*u*u/7.0
ELSE
  N = log(chi)/(2.*f)
ENDIF

UMean = (U_L+U_R)/(2.*N)
END SUBROUTINE getLogMean

!==================================================================================================================================
!> Computes the kinetic energy preservign surface flux
!==================================================================================================================================
PPURE SUBROUTINE KepALESurfaceFlux(U_LL,U_RR,MeshVel,F)
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_2Var),INTENT(IN)    :: U_LL,U_RR !< variables at the left-/right-Surfaces
REAL,                   INTENT(IN)    :: MeshVel   !< mesh velocity
REAL,DIMENSION(PP_nVar),INTENT(INOUT) :: F         !< flux, will be augmented by KEP ALE flux
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Kinetic energy stable ALE fluxes
F(1)=F(1)-0.5* MeshVel*(U_LL(DENS)+U_RR(DENS))
F(2)=F(2)-0.25*MeshVel*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL1)+U_RR(VEL1))
F(3)=F(3)-0.25*MeshVel*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL2)+U_RR(VEL2))
#if PP_dim == 3
F(4)=F(4)-0.25*MeshVel*(U_LL(DENS)+U_RR(DENS))*(U_LL(VEL3)+U_RR(VEL3))
#endif
F(5)=F(5)-0.25*MeshVel*(U_LL(DENS)+U_RR(DENS))*(U_LL(ENER)/U_LL(DENS)+U_RR(ENER)/U_RR(DENS))

END SUBROUTINE KepALESurfaceFlux

!==================================================================================================================================
!> Computes the kinetic energy preserving Split-Flux
!> Attention: Factor 2 from differentiation matrix is already been considered
!==================================================================================================================================
PPURE SUBROUTINE KepALEVolumeFlux(URef,UPrimRef,MeshVelRef,U,UPrim,MeshVel,fTilde,gTilde&
#if PP_dim == 3
                                  ,hTilde&
#endif
                                 )
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar    ),INTENT(IN)    :: URef,U             !< conserved variables
REAL,DIMENSION(PP_nVarPrim),INTENT(IN)    :: UPrimRef,UPrim     !< primitive variables
REAL,DIMENSION(1:3        ),INTENT(IN)    :: MeshVelRef,MeshVel !< mesh velocities
REAL,DIMENSION(PP_nVar)    ,INTENT(INOUT) :: fTilde,gTilde      !< flux in physical space
#if PP_dim == 3
REAL,DIMENSION(PP_nVar)    ,INTENT(INOUT) :: hTilde             !< flux in physical space
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
! Kinetic energy stable ALE fluxes
fTilde(1) = fTilde(1)-0.5*( MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))
fTilde(2) = fTilde(2)-0.25*(MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))
fTilde(3) = fTilde(3)-0.25*(MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))
#if PP_dim == 3
fTilde(4) = fTilde(4)-0.25*(MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))
#endif
fTilde(5) = fTilde(5)-0.25*(MeshVelRef(1)+MeshVel(1))*(URef(1)+U(1))*(URef(5)/URef(1)+U(5)/U(1))

gTilde(1) = gTilde(1)-0.5*( MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))
gTilde(2) = gTilde(2)-0.25*(MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))
gTilde(3) = gTilde(3)-0.25*(MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))
#if PP_dim == 3
gTilde(4) = gTilde(4)-0.25*(MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))
#endif
gTilde(5) = gTilde(5)-0.25*(MeshVelRef(2)+MeshVel(2))*(URef(1)+U(1))*(URef(5)/URef(1)+U(5)/U(1))

#if PP_dim == 3
hTilde(1) = hTilde(1)-0.5*( MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))
hTilde(2) = hTilde(2)-0.25*(MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))*(UPrimRef(2)+UPrim(2))
hTilde(3) = hTilde(3)-0.25*(MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))*(UPrimRef(3)+UPrim(3))
hTilde(4) = hTilde(4)-0.25*(MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))*(UPrimRef(4)+UPrim(4))
hTilde(5) = hTilde(5)-0.25*(MeshVelRef(3)+MeshVel(3))*(URef(1)+U(1))*(URef(5)/URef(1)+U(5)/U(1))
#endif

END SUBROUTINE KepALEVolumeFlux

END MODULE MOD_SplitFlux
