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
#include "flexi.h"

!===================================================================================================================================
!> Tool used to pre-calculate the distance to the closest solid wall for each volume gauss point when needed e.g. for the 
!> RANS SA turbulence model.
!> General process is as follows:
!>   * Read in of mesh (global, only single execution)
!>   * For each volume gauss points: 
!>       * Coarse search. Super-sampling of all solid wall surfaces to find the closest surface.
!>       * Minimization of the square of the distance function using a simple conjugate gradient algorithm on that specific surface.
!>         The coarse search is used as a starting point.
!===================================================================================================================================
PROGRAM SortIcingSides
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
#if USE_MPI
USE MOD_MPI,                     ONLY: FinalizeMPI
#endif
USE MOD_Mesh,                    ONLY: DefineParametersMesh,InitMesh,FinalizeMesh
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_SortIcingSides,          ONLY: CalcSortIcingSides
USE MOD_Interpolation,           ONLY: InitInterpolation,FinalizeInterpolation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=255) :: IniFile_dummy = ".dummy.ini"
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

! Define parameters needed
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMesh()

! check if mesh file is given
IF ((nArgs.NE.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'h5')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: SortIcingSides mesh.h5')
END IF
OPEN(1, FILE=IniFile_dummy)
CLOSE(1)
CALL prms%read_options(IniFile_dummy)

IF(PP_dim.EQ.2) CALL Abort(__STAMP__,'PLEASE COMPILE 3D')

SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') "     SORT ICING SIDES      "
SWRITE(UNIT_stdOut,'(132("="))')

! Initialization
CALL InitIOHDF5()
CALL InitInterpolation(NIn=1) !dummy
CALL InitMesh(meshMode=1,MeshFile_IN=Args(1),doDeallocateNodeCoords=.FALSE.)

CALL CalcSortIcingSides()

CALL FinalizeInterpolation()

CALL FinalizeMesh()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' SORT ICING SIDES TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM SortIcingSides
