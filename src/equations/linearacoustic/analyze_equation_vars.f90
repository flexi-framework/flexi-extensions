!=================================================================================================================================
! Copyright (c) 2010-2024  Prof. Claus-Dieter Munz
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
!> Contains global variables used by the equation specific analyze modules.
!==================================================================================================================================
MODULE MOD_AnalyzeEquation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL              :: doCalcTimeAverage   =.FALSE.      !< marks if time averaging should be performed
LOGICAL              :: doCalcAcousticEnergy=.FALSE.      !< marks if the mean acoustic energy should be computed
LOGICAL              :: doWriteAcousticEnergy=.FALSE.     !< marks if the mean acoustic energy should be written to file
CHARACTER(LEN=255)   :: FileName_AcousticEnergy           !< File Name for mean acoustic energy
!==================================================================================================================================
END MODULE MOD_AnalyzeEquation_Vars
