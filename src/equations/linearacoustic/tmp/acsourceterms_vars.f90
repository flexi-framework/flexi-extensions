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
!> Contains global variables used by the equation specific analyze modules.
!==================================================================================================================================
MODULE MOD_ACSourceterms_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
! Acoustic source terms

INTEGER              :: nCalcAcSources=1                  !< output acoustic sources every nth time step
LOGICAL              :: doCalcAcSources =.FALSE.          !< marks if acoustic source terms should be computed
REAL   ,ALLOCATABLE  :: AcSource(:,:,:,:,:)               !< time averaged solution U
LOGICAL,ALLOCATABLE  :: CalcAc(:)                         !< variables for which time averages should be computed (global indexing)
INTEGER,ALLOCATABLE  :: iAc(:)                            !< map from (global) VariableList to index in UAvg array
INTEGER              :: nVarAc                            !< number of time averag variables
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAcOut(:)        !< time averaged variable names
!==================================================================================================================================
END MODULE MOD_ACSourceterms_Vars
