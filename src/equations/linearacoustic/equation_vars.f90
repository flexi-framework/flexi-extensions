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
!> Contains the constant Advection Velocity Vector used for the linearized euler equations
!==================================================================================================================================
MODULE MOD_Equation_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: doCalcSource             !< automatically set by calcsource itself
INTEGER                  :: IniExactFunc
INTEGER                  :: IniRefState              !< RefState for initialization (case IniExactFunc=1 only)
INTEGER                  :: nRefState                !< number of refstates defined in parameter file
REAL,ALLOCATABLE         :: RefStateCons(:,:)        !< solution refstates (fluctuations to primitive state (rho,u,v,w,p). 
                                                     !< Cons only for compatibility reasons.
CHARACTER(LEN=255)       :: BCStateFile              !< file containing the reference solution on the boundary to be used as BC

! Boundary condition arrays
REAL,ALLOCATABLE         :: BCData(:,:,:,:)          !< Buffer array for BC data
INTEGER,ALLOCATABLE      :: nBCByType(:)             !< Number of sides for each boundary
INTEGER,ALLOCATABLE      :: BCSideID(:,:)            !< SideIDs for BC types

REAL                     :: Kappa                    !< heat capacity ratio / isentropic exponent
REAL                     :: KappaM1                  !< = $\kappa - 1$
REAL                     :: sKappaM1                 !< = $1/(\kappa -1)$
REAL                     :: KappaP1                  !< = $\kappa + 1$
REAL                     :: sKappaP1                 !< = $1/(\kappa +1)$
REAL                     :: s43,s23

! Acoustic source term
LOGICAL                  :: useAcSources
LOGICAL                  :: readAcSourcesFromFile
REAL,ALLOCATABLE,TARGET  :: source(:,:,:,:,:)

CHARACTER(LEN=255),DIMENSION(5),PARAMETER :: StrVarNames =&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure'/) !< solution variable names
CHARACTER(LEN=255),DIMENSION(5),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure'/) !< primitive variable names

CHARACTER(LEN=255),DIMENSION(6),PARAMETER :: StrVarNamesBase =&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Soundspeed'/) !< solution variable names

LOGICAL           :: EquationInitIsDone=.FALSE.!< Init switch  
!==================================================================================================================================
END MODULE MOD_Equation_Vars
