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
!> Contains the surface integral for the GCL. Using the poor man's template, the routine included for the main equations is reused.
!==================================================================================================================================
MODULE MOD_GCL_SurfInt
IMPLICIT NONE
PRIVATE

#undef WITHnVars
INTEGER,PARAMETER :: TP_nVar = 1 !< For the GCL, we consider a scalar quantity

INTERFACE GCL_SurfInt
  MODULE PROCEDURE SurfInt
END INTERFACE

INTERFACE GCL_DoSurfInt
  MODULE PROCEDURE DoSurfInt
END INTERFACE

PUBLIC::GCL_SurfInt,GCL_DoSurfInt

CONTAINS
#include "../dg/surfint.t90"
END MODULE MOD_GCL_SurfInt
