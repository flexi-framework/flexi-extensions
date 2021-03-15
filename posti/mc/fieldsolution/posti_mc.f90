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


PROGRAM MC
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_MC_Vars
USE MOD_MC                    ,ONLY: DefineParametersMC, InitMC,  ComputeMCEstimator, FinalizeMC
USE MOD_MC_Output             ,ONLY: WriteMeanAndVarianceMCToHDF5
USE MOD_ReadInTools
USE MOD_StringTools             ,ONLY: STRICMP,GetFileExtension
USE MOD_MPI                     ,ONLY: DefineParametersMPI,InitMPI
USE MOD_IO_HDF5                 ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
#if USE_MPI
USE MOD_MPI                     ,ONLY: InitMPIvars,FinalizeMPI
#endif
USE MOD_Eos,               ONLY:DefineParametersEos,InitEOS
!USE MOD_Equation,          ONLY:DefineParametersEquation,InitEquation
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
CALL SetStackSizeUnlimited()
CALL InitMPI()
IF (nProcessors.GT.1) CALL CollectiveStop(__STAMP__, &
     'This tool is designed only for single execution!')

CALL ParseCommandlineArguments()

WRITE(UNIT_stdOut,'(A)') " ||===========================================================================================================||"
WRITE(UNIT_stdOut,'(A)') " || Compute MC estimate                                                                                     ! ||"
WRITE(UNIT_stdOut,'(A)') " ||===========================================================================================================||"
WRITE(UNIT_stdOut,'(A)')

! Define parameters needed
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersMC()
CALL DefineParametersEos()
!CALL DefineParametersEquation()
CALL prms%read_options(Args(1))

IF (doPrintHelp.GT.0) THEN
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.1).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: postiMC prm-file statefiles')
END IF
CALL InitIOHDF5()
CALL InitEOS()
!CALL InitEquation()
CALL InitMC()
nStateFiles=nArgs-1
!----------------------------------------------------------------------------------------------------------------------------------!
! Compute Modes, Mean and Variance
!----------------------------------------------------------------------------------------------------------------------------------!
CALL ComputeMCEstimator()

SWRITE(UNIT_stdOut,'(A)') ' WRITING  MEAN and VARIANCE...'
CALL WriteMeanAndVarianceMCToHDF5()
CALL FinalizeMC()
#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
CALL FinalizeMPI()
#endif
WRITE(UNIT_stdOut,'(132("="))')
WRITE(UNIT_stdOut,'(A)') ' NISP TOOL FINISHED! '
WRITE(UNIT_stdOut,'(132("="))')

END PROGRAM MC
