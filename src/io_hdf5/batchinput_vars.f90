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
!> Contains variables relevant for indicators
!==================================================================================================================================
MODULE MOD_BatchInput_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)             :: StochFile

INTEGER                        :: nStochVars
REAL,ALLOCATABLE               :: iStochSample(:)
CHARACTER(LEN=255),ALLOCATABLE :: StochVarNames(:)
INTEGER,ALLOCATABLE            :: iOccurrence(:)
INTEGER,ALLOCATABLE            :: iArray(:)

INTEGER                        :: nLevelVarsInt
INTEGER                        :: nLevelVarsReal
INTEGER                        :: nLevelVarsStr
CHARACTER(LEN=255),ALLOCATABLE :: LevelVarNamesInt(:)
CHARACTER(LEN=255),ALLOCATABLE :: LevelVarNamesReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: LevelVarNamesStr(:)
INTEGER,ALLOCATABLE            :: LevelVarsInt(:)
REAL,ALLOCATABLE               :: LevelVarsReal(:)
CHARACTER(LEN=255),ALLOCATABLE :: LevelVarsStr(:)

LOGICAL                        :: BatchMode
LOGICAL                        :: BatchInputInitIsDone=.FALSE.
LOGICAL                        :: StochMeshFileExists=.FALSE.
CHARACTER(LEN=255)             :: StochMeshFile
!==================================================================================================================================
END MODULE MOD_BatchInput_Vars
