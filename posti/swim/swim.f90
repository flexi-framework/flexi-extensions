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
!>
!===================================================================================================================================
PROGRAM Swim
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Swim_Vars
USE MOD_Commandline_Arguments
USE MOD_ReadInTools
USE MOD_StringTools,             ONLY: STRICMP,GetFileExtension
USE MOD_MPI,                     ONLY: DefineParametersMPI,InitMPI
USE MOD_IO_HDF5,                 ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI,                     ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_SwimSolver,              ONLY: InitSwim
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

! Parse parameters
! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.EQ.1).AND.STRICMP(GetFileExtension(Args(1)),'h5')) THEN
  dointerpolatex = .FALSE.
  StateFile = Args(1)
ELSEIF((nArgs.EQ.5).AND.STRICMP(Args(1),"--InterpolateX").AND.STRICMP(GetFileExtension(Args(5)),'h5')) THEN
  dointerpolatex = .TRUE.
  READ(Args(2),*) xMin
  READ(Args(3),*) xMax
  READ(Args(4),*) nPts
  StateFile = Args(5)
ELSE
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: swim statefile')
END IF
OPEN(1, FILE=IniFile_dummy)
CLOSE(1)
CALL prms%read_options(IniFile_dummy)


SWRITE(UNIT_stdOut,'(132("="))')
SWRITE(UNIT_stdOut,'(A)') "     SWIM     "
SWRITE(UNIT_stdOut,'(132("="))')

! Initialization
CALL InitIOHDF5()

CALL InitSwim()

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif

WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' SWIM FINISHED!'
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM Swim
