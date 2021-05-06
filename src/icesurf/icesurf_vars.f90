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
!==================================================================================================================================
!> Contains global variables used by the DG modules.
!==================================================================================================================================
MODULE MOD_IceSurf_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL             :: doCalcIceSurfData = .FALSE.

LOGICAL             :: doAvgIceSurf = .TRUE.
LOGICAL             :: SurfFirstZOnly = .TRUE.
REAL,ALLOCATABLE    :: IceSurfData(:,:,:,:)
REAL,ALLOCATABLE    :: UAvg_master(:,:,:,:)
REAL,ALLOCATABLE    :: UAvg_slave(:,:,:,:)
INTEGER,ALLOCATABLE :: MapIceSurf(:)

INTEGER             :: nWallSides
INTEGER             :: nWallSidesglob
INTEGER             :: offsetWallSides
INTEGER             :: NOutSurf
INTEGER             :: nVarSurf

! Orientation of P and Q relative to TangVec1 and TangVec2 (number => which vector, sign => which orientation)
! Derivation below
#if PP_dim == 3
INTEGER             :: PDir_T(6) =(/-2,-2, 1, 2,-2, 1/)
INTEGER             :: QDir_T(6) =(/ 1, 1, 2, 1, 1, 2/)
#else
INTEGER             :: PDir_T(6) =(/ 0, 1, 1,-1,-1, 0/)
#endif
! Derivation (3D): 
! Direction and orientation of p and q (from CGNS cube): 
! PDir =  (/ ETA_PLUS,   XI_PLUS,    ETA_PLUS,  XI_MINUS,  ZETA_PLUS,  XI_PLUS   /)
! QDir =  (/ XI_PLUS,    ZETA_PLUS,  ZETA_PLUS, ZETA_PLUS, ETA_PLUS,   ETA_PLUS  /) 
! Direction and orientation of vectors: 
! - normal is just iLocSide (for order, see flexi.h) 
! - T1 from TangDirs in mesh_vars-f90 (always positive)
! - T2 from right-hand system / cross product T2 = N x T1 (use fingers with CGNS cube)
! NDir =  (/ ZETA_MINUS, ETA_MINUS,  XI_PLUS,   ETA_PLUS,  XI_MINUS,   ZETA_PLUS /)
! T1Dir = (/ XI_PLUS,    ZETA_PLUS,  ETA_PLUS,  ZETA_PLUS, ETA_PLUS,   XI_PLUS   /)
! T2Dir = (/ ETA_MINUS,  XI_MINUS,   ZETA_PLUS, XI_MINUS,  ZETA_MINUS, ETA_PLUS  /)
! Comparison yields PDir_T and QDir_T

! Derivation (2D): 
! - NDir is as in 3D
! - PDir is always anti-clockwise (see for example CGNS_VolToSide in mappings.f90)
! - T1 is always positive
! - T2 is always zero and q = 0, only orientation  p += T1 matters 
! - comparinson yields PDir_T
! NDir =  (/ 0, ETA_MINUS, XI_PLUS,  ETA_PLUS, XI_MINUS,  0 /)
! PDir =  (/ 0, XI_PLUS,   ETA_PLUS, XI_MINUS, ETA_MINUS, 0 /)
! T1Dir = (/ 0, XI_PLUS,   ETA_PLUS, XI_PLUS,  ETA_PLUS,  0 /)

LOGICAL             :: IceSurfInitIsDone = .FALSE.
END MODULE MOD_IceSurf_Vars
