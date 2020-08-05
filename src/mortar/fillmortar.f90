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
!> \brief Auxiliary routines to interpolate the normal vectors of the sliding mesh mortars.
!> 
!> Contains the routines to
!> - interpolate the normal and tangential vectors at the sm interface onto the sm mortars
!> - flip the local orientation of the sides across the sliding mesh interface
!==================================================================================================================================
MODULE MOD_MortarSM
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = 3

INTERFACE DoFlip3
  MODULE PROCEDURE DoFlip
END INTERFACE

INTERFACE InterpolateSM3
  MODULE PROCEDURE InterpolateSM
END INTERFACE

INTERFACE ProjectSM3
  MODULE PROCEDURE ProjectSM
END INTERFACE

PUBLIC:: DoFlip3, InterpolateSM3, ProjectSM3

CONTAINS
#include "fillmortar_sm.t90"
END MODULE MOD_MortarSM

!==================================================================================================================================
!> \brief Routines that perform the projection operation between nonconforming interfaces using the operators set up in module
!> mortar
!>
!> Contains the routines to
!> - interpolate the solution at the large sides to the small ones, which are used for flux computation
!> - project the flux from the small sides back to the large ones
!==================================================================================================================================
MODULE MOD_FillMortar
IMPLICIT NONE
PRIVATE

#define WITHnVar 1

INTERFACE U_Mortar
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

INTERFACE U_MortarSM
  MODULE PROCEDURE U_MortarSM
END INTERFACE

INTERFACE Flux_MortarSM
  MODULE PROCEDURE Flux_MortarSM
END INTERFACE

PUBLIC::U_Mortar,Flux_Mortar,U_MortarSM,Flux_MortarSM

CONTAINS
#include "fillmortar.t90"
#include "fillmortar_sm.t90"
END MODULE MOD_FillMortar

!==================================================================================================================================
!> Routines that perform the projection operation between nonconforming interfaces of conservative variables
!==================================================================================================================================
MODULE MOD_FillMortarCons
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVar

INTERFACE U_MortarCons
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_MortarCons
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

INTERFACE U_MortarConsSM
  MODULE PROCEDURE U_MortarSM
END INTERFACE

INTERFACE Flux_MortarConsSM
  MODULE PROCEDURE Flux_MortarSM
END INTERFACE

PUBLIC::U_MortarCons,Flux_MortarCons,U_MortarConsSM,Flux_MortarConsSM

CONTAINS
#include "fillmortar.t90"
#include "fillmortar_sm.t90"
END MODULE MOD_FillMortarCons

!==================================================================================================================================
!> Routines that perform the projection operation between nonconforming interfaces of primitive variables
!==================================================================================================================================
MODULE MOD_FillMortarPrim
IMPLICIT NONE
PRIVATE

#undef WITHnVar
INTEGER,PARAMETER :: TP_nVar = PP_nVarPrim

INTERFACE U_MortarPrim
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_MortarPrim
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

INTERFACE U_MortarPrimSM
  MODULE PROCEDURE U_MortarSM
END INTERFACE

INTERFACE Flux_MortarPrimSM
  MODULE PROCEDURE Flux_MortarSM
END INTERFACE

PUBLIC::U_MortarPrim,Flux_MortarPrim,U_MortarPrimSM,Flux_MortarPrimSM

CONTAINS
#include "fillmortar.t90"
#include "fillmortar_sm.t90"
END MODULE MOD_FillMortarPrim

!==================================================================================================================================
!> Routines that perform the projection operation between nonconforming interfaces of a scalar variable
!==================================================================================================================================
MODULE MOD_FillMortar1
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: TP_nVar = 1

INTERFACE U_Mortar1
  MODULE PROCEDURE U_Mortar
END INTERFACE

INTERFACE Flux_Mortar1
  MODULE PROCEDURE Flux_Mortar
END INTERFACE

INTERFACE U_MortarSM1
  MODULE PROCEDURE U_MortarSM
END INTERFACE

INTERFACE Flux_MortarSM1
  MODULE PROCEDURE Flux_MortarSM
END INTERFACE

PUBLIC::U_Mortar1,Flux_Mortar1,U_MortarSM1,Flux_MortarSM1

CONTAINS
#include "fillmortar.t90"
#include "fillmortar_sm.t90"
END MODULE MOD_FillMortar1
