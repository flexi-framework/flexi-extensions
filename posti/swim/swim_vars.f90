!=================================================================================================================================
! Copyright (c) 2016  Prof. Claus-Dieter Munz 
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

!===================================================================================================================================
!> Contains global variables provided by the SortIcingSides routines
!===================================================================================================================================
MODULE MOD_Swim_Vars
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PUBLIC
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
CHARACTER(LEN=255)             :: MeshFile      !< name of meshfile (contains indexing)
CHARACTER(LEN=255)             :: StateFile     !< name of hdf5 statefile (contains surface data)
LOGICAL             :: dummy
INTEGER             :: nWallSides
INTEGER             :: nCells
INTEGER,ALLOCATABLE :: Mapping(:,:)
REAL,ALLOCATABLE    :: VolData(:,:,:,:,:)
INTEGER             :: ICS_N,ICS_NZ
INTEGER             :: WriteDim
INTEGER             :: nRuns
INTEGER             :: nVarSurf
REAL,ALLOCATABLE    :: XY(:,:)
REAL,ALLOCATABLE    :: SwimData(:,:,:)
!===================================================================================================================================
END MODULE MOD_Swim_Vars
