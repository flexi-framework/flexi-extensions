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
#include "flexi.h"

MODULE MOD_VisualizeRP
IMPLICIT NONE
PRIVATE

INTERFACE DefineParameters
  MODULE PROCEDURE DefineParameters
END INTERFACE
  
INTERFACE InitParameters
  MODULE PROCEDURE InitParameters
END INTERFACE

INTERFACE Build_mapCalc_mapVisu
  MODULE PROCEDURE Build_mapCalc_mapVisu
END INTERFACE

PUBLIC::DefineParameters,InitParameters,Build_mapCalc_mapVisu
CONTAINS

!===================================================================================================================================
!> Initialize parameter variables of Posti tool
!===================================================================================================================================
SUBROUTINE DefineParameters()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!===================================================================================================================================
CALL prms%SetSection('Visualize Record Points')

CALL prms%CreateStringOption( 'ProjectName'        ,"Name of the project")
CALL prms%CreateStringOption( 'GroupName'          ,"Name(s) of the group(s) to visualize, must be equal to the name given in&
                                                     & preparerecordpoints tool",multiple=.TRUE.)
CALL prms%CreateStringOption( 'VarName'            ,"Variable name to visualize",multiple=.TRUE.)
CALL prms%CreateStringOption( 'RP_DefFile'         ,"Path to the *RPset.h5 file")

CALL prms%CreateLogicalOption('usePrims'           ,"Set to indicate that the RP file contains the primitive and not the&
                                                    & conservative variables",".FALSE.")

CALL prms%CreateLogicalOption('OutputTimeData'     ,"Should the time series be written? Not compatible with TimeAvg and FFT&
                                                     & options!",".FALSE.")
CALL prms%CreateLogicalOption('OutputTimeAverage'  ,"Should the time average be computed and written?",".FALSE.")
CALL prms%CreateLogicalOption('doFluctuations'     ,"Should the fluctuations be computed and written?",".FALSE.")
CALL prms%CreateLogicalOption('equiTimeSpacing'    ,"Set to interpolate the temporal data to equdistant time steps&
                                                     & (always done for operations requiring FFTs)",".FALSE.")

CALL prms%CreateLogicalOption('OutputPoints'       ,"General option to turn off the output of points",".TRUE.")
CALL prms%CreateLogicalOption('OutputLines'        ,"General option to turn off the output of lines",".TRUE.")
CALL prms%CreateLogicalOption('OutputPlanes'       ,"General option to turn off the output of planes",".TRUE.")

CALL prms%CreateLogicalOption('doFFT'              ,"Calculate a fast Fourier transform of the time signal",".FALSE.")
CALL prms%CreateLogicalOption('doPSD'              ,"Calculate the power spectral density of the time signal",".FALSE.")
CALL prms%CreateIntOption    ('nBlocks'            ,"Specify the number of blocks over the time signal used for spectral averaging&
                                                     & when calculating spectral quantities",'1')
CALL prms%CreateRealOption   ('SamplingFreq'       ,"Instead of specifying the number of blocks, the sampling frequency in&
                                                     & combination with the block size can be set - the number of blocks&
                                                     & will then be calculated.")
CALL prms%CreateIntOption    ('BlockSize'          ,"Size of the blocks (in samples) if sampling frequency is given")
CALL prms%CreateRealOption   ('CutoffFreq'         ,"Specify smallest considered frequency in spectral analysis")
CALL prms%CreateLogicalOption('hanning'            ,"Set to use the Hann window when performing spectral analysis",".FALSE.")
CALL prms%CreateLogicalOption('fourthDeriv'        ,"Set to calculate the fourth derivative of the time signal (signal will be&
                                                     & interpolated to a coarse grid since truncation error is large",".FALSE.")
CALL prms%CreateLogicalOption('ThirdOct'           ,"TODO",".FALSE.")
CALL prms%CreateRealOption   ('u_inf'              ,"TODO")
CALL prms%CreateRealOption   ('chord'              ,"TODO")

CALL prms%CreateLogicalOption('doTurb'             ,"Set to compute a temporal FFT for each RP and compute turbulent quantities&
                                                    & like the kinetic energy over wave number",".FALSE.")

CALL prms%CreateLogicalOption('Plane_doBLProps'    ,"Set to calculate seperate boundary layer quantities for boundary layer&
                                                     & planes",".FALSE.")
CALL prms%CreateIntOption    ('Plane_BLvelScaling' ,"Choose scaling for boundary layer quantities. 0: no scaling, 1: laminar&
                                                     & scaling, 3: turbulent scaling")
CALL prms%CreateLogicalOption('Plane_LocalCoords'  ,"Set to use local instead of global coordinates along planes",".FALSE.")
CALL prms%CreateLogicalOption('Plane_LocalVel'     ,"Set to use local instead of global velocities along planes",".FALSE.")

CALL prms%CreateLogicalOption('Line_LocalCoords'   ,"Set to use local instead of global coordinates along lines",".FALSE.")
CALL prms%CreateLogicalOption('Line_LocalVel'      ,"Set to use local instead of global velocities along lines",".FALSE.")
CALL prms%CreateRealArrayOption('Line_LocalVel_vec',"Vector used for local velocity computation along line")

CALL prms%CreateLogicalOption('doFilter'           ,"Set to perform temporal filtering for each RP",".FALSE.")
CALL prms%CreateRealOption   ('FilterWidth'        ,"Width of the temporal filter")
CALL prms%CreateIntOption    ('FilterMode'         ,"Set to 0 for low pass filter and to 1 for high pass filter")

CALL prms%CreateRealOption   ('mu0'                ,"Kinematic viscosity, needed for turbulent quantities")

CALL prms%CreateIntOption    ('SkipSample'         ,"Used to skip every n-th RP evaluation")
CALL prms%CreateIntOption    ('OutputFormat'       ,"Choose the main format for output. 0: ParaView, 2: HDF5")
END SUBROUTINE DefineParameters

!===================================================================================================================================
!> Read parameters of Posti tool
!===================================================================================================================================
SUBROUTINE InitParameters()
! MODULES
USE MOD_Globals
USE MOD_Readintools         ,ONLY:GETINT,GETREAL,GETLOGICAL,GETSTR,GETREALARRAY,CountOption
USE MOD_ParametersVisu
USE MOD_RPInterpolation_Vars,ONLY:calcTimeAverage
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                :: iGroup
!===================================================================================================================================
Projectname=GETSTR('ProjectName')

! =============================================================================== !
! RP INFO
! =============================================================================== !
nGroups_visu=CountOption('GroupName')
ALLOCATE(GroupNames_visu(nGroups_visu))
DO iGroup=1,nGroups_visu
  GroupNames_visu(iGroup)=GETSTR('Groupname','none')
END DO
RP_SET_defined=.FALSE.
RP_DefFile=GETSTR('RP_DefFile','none')
IF(TRIM(RP_defFile).NE.'none') THEN
  RP_SET_defined=.TRUE.
END IF

! use primitive variables for derived quantities if they exist in the state file
usePrims=GETLOGICAL('usePrims','.FALSE.')

! =============================================================================== !
! TIME INTERVAL
! =============================================================================== !
OutputTimeData   =GETLOGICAL('OutputTimeData','.FALSE.')
OutputTimeAverage=GETLOGICAL('OutputTimeAverage','.FALSE.')
doFluctuations   =GETLOGICAL('doFluctuations','.FALSE.')
IF(OutputTimeAverage .OR. doFluctuations) CalcTimeAverage=.TRUE.
IF(doFluctuations) OutputTimeData = .TRUE.

doFilter=GETLOGICAL('doFilter','.FALSE.')
IF(doFilter) THEN
  FilterWidth=GETREAL('FilterWidth')
  FilterMode=GETINT('FilterMode','0')
END IF
! =============================================================================== !
! FOURIER TRANSFORM
! =============================================================================== !
doFFT=GETLOGICAL('doFFT','.FALSE.')
doPSD=GETLOGICAL('doPSD','.FALSE.')
IF(doFFT.OR.doPSD) doSpec=.TRUE.
IF(doSpec) THEN
  ! two readin "modes" for spectrum averaging:
  ! 1. Prescription of number of blocks
  nBlocks=GETINT('nBlocks','1')
  ! 2. Prescription of Sampling Frequency and Blocksize
  samplingFreq=GETREAL('SamplingFreq','-999')
  IF(samplingFreq.GT.0.) THEN
    BlockSize=GETINT('BlockSize')
  END IF
  doHanning=GETLOGICAL('hanning','.FALSE.')
  fourthDeriv=GETLOGICAL('fourthDeriv','.FALSE.')
  ThirdOct=GETLOGICAL('ThirdOct','.FALSE.')
  IF (ThirdOct) THEN
    u_infPhys   = GETREAL('u_inf') !velocity for re-dimensionalization of frequency
    chordPhys = GETREAL('chord')   !length for re-dimensionalization of frequency
  END IF
END IF

doTurb=GETLOGICAL('doTurb','.FALSE.')

IF(doSpec .OR. doTurb) cutoffFreq=GETREAL('cutoffFreq','-999.')

! for any FFT operation we need equidistant time spacing
equiTimeSpacing=.FALSE.
IF(OutputTimeData) equiTimeSpacing=GETLOGICAL('equiTimeSpacing','.FALSE.')
IF(doTurb.OR.doSpec) equiTimeSpacing=.TRUE.

! =============================================================================== !
! PLANE OPTIONS
! =============================================================================== !
OutputPlanes      =GETLOGICAL('OutputPlanes','.TRUE.')
Plane_LocalCoords =GETLOGICAL('Plane_LocalCoords','.FALSE.')
Plane_LocalVel    =GETLOGICAL('Plane_LocalVel','.FALSE.')
Plane_doBLProps   =GETLOGICAL('Plane_doBLProps','.FALSE.')
IF(Plane_doBLProps) THEN ! for BL properties we need local coords and velocities
  WRITE(UNIT_StdOut,'(A)')' BL properties depend on local velocities and coordinates'
  WRITE(UNIT_StdOut,'(A)')' and are calculated based on time-averaged data.'
  WRITE(UNIT_StdOut,'(A)')' Setting Plane_localCoords=.TRUE. and Plane_localVel=.TRUE..'
  CalcTimeAverage  =.TRUE.
  OutputTimeAverage=.TRUE.
  Plane_LocalCoords=.TRUE.
  Plane_LocalVel   =.TRUE.
  Plane_BLvelScaling  =GETINT('Plane_BLvelScaling','0') ! 0 - no scaling.
  ! 1 - "laminar scaling": scale velocity with u_delta and PlaneY with delta99
  ! 2 - "turbulent scaling:" calculate u+ and y+
END IF
! =============================================================================== !
! LINE OPTIONS
! =============================================================================== !
OutputLines      =GETLOGICAL('OutputLines',     '.TRUE.')
Line_LocalCoords =GETLOGICAL('Line_LocalCoords','.FALSE.')
Line_LocalVel    =GETLOGICAL('Line_LocalVel',   '.FALSE.')
IF(Line_LocalVel) THEN
  Line_LocalVel_vec=GETREALARRAY('Line_LocalVel_vec',3)
  ! normalize the vector
  Line_LocalVel_vec = Line_LocalVel_vec/NORM2(Line_LocalVel_vec)
END IF
! =============================================================================== !
! POINT OPTIONS
! =============================================================================== !
OutputPoints     =GETLOGICAL('OutputPoints','.TRUE.')

skip = GETINT('SkipSample','1')

IF(doTurb.OR.Plane_doBLProps) Mu0=GETREAL('Mu0')

OutputFormat = GETINT('OutputFormat','0')

! when no time averaging and no fft shall be performed, set the output signal to true
IF(.NOT.ANY( (/doSpec, OutputTimeAverage , doTurb /)) ) OutputTimeData=.TRUE.
END SUBROUTINE InitParameters


!===================================================================================================================================
!> This routine builds the mappings from the total number of variables available for visualization to number of calculation
!> and visualization variables.
!>  1. Read 'VarName' options from the parameter file. This are the quantities that will be visualized.
!>  2. Initialize the dependecy table
!>  3. check wether gradients are needed for any quantity. If this is the case, remove the conservative quantities from the
!>     dependecies of the primitive quantities (the primitive quantities are available directly, since the DGTimeDerivative_weakForm
!>     will be executed.
!>  4. build the 'mapCalc' that holds for each quantity that will be calculated the index in 'UCalc' array (0 if not calculated)
!>  5. build the 'mapVisu' that holds for each quantity that will be visualized the index in 'UVisu' array (0 if not visualized)
!===================================================================================================================================
SUBROUTINE Build_mapCalc_mapVisu()
USE MOD_ParametersVisu
USE MOD_ReadInTools     ,ONLY: GETSTR,CountOption
USE MOD_StringTools     ,ONLY: STRICMP
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLESIABLES
INTEGER             :: iVar,iVar2
CHARACTER(LEN=20)   :: format
!===================================================================================================================================
! Read Varnames from parameter file and fill
!   mapVisu = map, which stores at position x the position/index of the x.th quantity in the UVisu array
!             if a quantity is not visualized it is zero
ALLOCATE(mapVisu(1:nVarDep))
mapVisu = 0
! Compare varnames that should be visualized with availabe varnames
DO iVar=1,nVarVisu
  DO iVar2=1,nVarDep
    IF (STRICMP(VarNameVisu(iVar), VarNamesAll(iVar2))) THEN
      mapVisu(iVar2) = iVar
    END IF
  END DO
END DO

! Calculate all dependencies:
! For each quantity copy from all quantities that this quantity depends on the dependencies.
DO iVar=1,nVarDep
  DepTable(iVar,iVar) = 1
  DO iVar2=1,iVar-1
    IF (DepTable(iVar,iVar2).EQ.1) &
      DepTable(iVar,:) = MAX(DepTable(iVar,:), DepTable(iVar2,:))
  END DO
END DO

! Build :
!   mapCalc = map, which stores at position x the position/index of the x.th quantity in the UCalc array
!             if a quantity is not calculated it is zero
ALLOCATE(mapCalc(1:nVarDep))
mapCalc = 0
DO iVar=1,nVarDep
  IF (mapVisu(iVar).GT.0) THEN
    mapCalc = MAX(mapCalc,DepTable(iVar,1:nVarDep))
  END IF
END DO
! enumerate mapCalc
nVarCalc = 0
DO iVar=1,nVarDep
  IF (mapCalc(iVar).GT.0) THEN
    nVarCalc = nVarCalc + 1
    mapCalc(iVar) = nVarCalc
  END IF
END DO

! print the dependecy table
WRITE(format,'(I2)') SIZE(DepTable,2)
DO iVar=1,nVarDep
  WRITE (*,'('//format//'I2,A)') DepTable(iVar,:), " "//TRIM(VarNamesAll(iVar))
END DO

! print the mappings
WRITE(format,'(I2)') nVarDep
WRITE (*,'(A,'//format//'I3)') "mapCalc ",mapCalc
WRITE (*,'(A,'//format//'I3)') "mapVisu ",mapVisu

END SUBROUTINE Build_mapCalc_mapVisu

END MODULE MOD_VisualizeRP
