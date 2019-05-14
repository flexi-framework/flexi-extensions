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
! Add comments please!
!===================================================================================================================================


PROGRAM NISP
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Nisp_Vars
USE MOD_Nisp                    ,ONLY: DefineParametersNisp, InitNisp, ComputeModes, FinalizeNisp
USE MOD_Nisp_Output             ,ONLY: WriteMeanAndVarianceToHDF5
USE MOD_ReadInTools
USE MOD_IO_HDF5                 ,ONLY: DefineParametersIO_HDF5
USE MOD_StringTools             ,ONLY: STRICMP,GetFileExtension
!#ifdef MPI
USE MOD_MPI                     ,ONLY: DefineParametersMPI,InitMPI
!#endif /* MPI */
USE MOD_EOS_Vars                ,ONLY: KappaM1,R, Kappa
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI() ! NO PARALLELIZATION, ONLY FOR COMPILING WITH MPI FLAGS ON SOME MACHINES OR USING MPI-DEPENDANT HDF5
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)') " || Compute NISP modes                 ! ||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

! Define parameters needed
CALL DefineParametersMPI()
CALL DefineParametersNisp()


IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: postiNisp prm-file StochInput.h5 statefile')
END IF
CALL InitNisp()
!----------------------------------------------------------------------------------------------------------------------------------!
! Compute Modes, Mean and Variance
!----------------------------------------------------------------------------------------------------------------------------------!
CALL ComputeModes()
SWRITE(UNIT_stdOut,'(A)') ' WRITING  MEAN and VARIANCE...'
CALL WriteMeanAndVarianceToHDF5()
CALL FinalizeNisp()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' NISP TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM NISP
