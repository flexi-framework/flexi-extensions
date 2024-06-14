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

!==================================================================================================================================
!> Contains the global variables needed for the exact funtctions of the LinAdv equation system.
!==================================================================================================================================
MODULE MOD_Exactfunc_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!----------------------------------------------------------------------------------------------------------------------------------
REAL              :: xP(3)                      !< Pressure pulse/source initial position CASE(3,4,5,6) (x,y,z)
REAL              :: sigmaSq                    !< Initial pulse/source width  CASE(3,4,5)
REAL              :: f                          !< frequency for monopole source CASE(3,6)
REAL              :: Amp                        !< amplitude for monopole source CASE(3,6)
!==================================================================================================================================

END MODULE MOD_Exactfunc_Vars
