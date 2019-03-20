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
!> Add comments please!
!===================================================================================================================================
PROGRAM Nisp_RP
! MODULES
USE MOD_Globals
USE MOD_Commandline_Arguments
USE MOD_Nisp_RP_Vars
USE MOD_Nisp_RP                     ,ONLY: DefineParametersNisp_RP,InitNisp_RP,PerformSampleFFT,ComputeModes,FinalizeNisp_RP
USE MOD_Nisp_RP_Output              ,ONLY: WriteMeanAndVarianceToHDF5!,SamplesToHDF5
!USE MOD_ReadInTools
USE MOD_IO_HDF5                     ,ONLY: DefineParametersIO_HDF5,InitIOHDF5
USE MOD_StringTools                 ,ONLY: STRICMP,GetFileExtension
USE MOD_ReadInTools
!#ifdef MPI
USE MOD_MPI                         ,ONLY: DefineParametersMPI,InitMPI
!#endif /* MPI */
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
WRITE(UNIT_stdOut,'(A)') " || Compute NISP modes based on RP data! ||"
WRITE(UNIT_stdOut,'(A)') " ||======================================||"
WRITE(UNIT_stdOut,'(A)')

! Define parameters needed
CALL DefineParametersMPI()
CALL DefineParametersIO_HDF5()
CALL DefineParametersNisp_RP()

! check for command line argument --help or --markdown
IF (doPrintHelp.GT.0) THEN
WRITE(UNIT_stdOut,'(A)') " || Compute NISP modes based on RP data!    ||"
  CALL PrintDefaultParameterFile(doPrintHelp.EQ.2, Args(1))
  STOP
END IF
! check if parameter file is given
IF ((nArgs.LT.3).OR.(.NOT.(STRICMP(GetFileExtension(Args(1)),'ini')))) THEN
  CALL CollectiveStop(__STAMP__,'ERROR - Invalid syntax. Please use: postiNisp_RP prm-file StochInput.h5 statefile [statefiles]')
END IF

CALL prms%read_options(Args(1))

!======================================================
! Init
!======================================================
CALL InitIOHDF5()
CALL InitNisp_RP()

!======================================================
! Main
!======================================================
CALL PerformSampleFFT()
CALL ComputeModes()
CALL WriteMeanAndVarianceToHDF5()
!CALL WriteSamplesToHDF5()

!======================================================
! Finalize
!======================================================
CALL FinalizeNisp_RP()
#ifdef MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) &
  CALL abort(__STAMP__,'MPI finalize error',iError)
#endif
END PROGRAM Nisp_RP
