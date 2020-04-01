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

!==================================================================================================================================
!> Control program of the Flexi code.
!==================================================================================================================================
PROGRAM FlexiBatch
! MODULES
USE MOD_Globals
USE MOD_Flexi
USE MOD_TimeDisc   ,ONLY: TimeDisc
USE MOD_IO_HDF5    ,ONLY: OpenDataFile,CloseDataFile,File_ID,InitMPIInfo
USE MOD_HDF5_Input ,ONLY: ReadAttribute,DatasetExists
USE MOD_BatchInput_Vars, ONLY: StochFile,BatchMode
USE MOD_StringTools ,ONLY: STRICMP,GetFileExtension
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgsLoc,iArg,Color,nArgsAdd
CHARACTER(LEN=255),ALLOCATABLE :: ArgsLoc(:)
CHARACTER(LEN=255)             :: FirstArg
LOGICAL                        :: exists,isActive
REAL                           :: GlobalStartTime,GlobalEndTime
!==================================================================================================================================
#if USE_MPI
CALL MPI_INIT(iError)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myGlobalRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nGlobalProcessors, iError)

MPI_COMM_ACTIVE=MPI_COMM_WORLD
#else
myGlobalRank=0
nGlobalProcessors=1
#endif
MPIGlobalRoot=(myGlobalRank .EQ. 0)
CALL InitMPIInfo()


! read command line arguments
nArgsLoc = COMMAND_ARGUMENT_COUNT()
IF(nArgsLoc.LT.1) CALL Abort(__STAMP__,'Provide Ini File!')
CALL GET_COMMAND_ARGUMENT(1,FirstArg)
BatchMode = STRICMP(GetFileExtension(FirstArg),'h5')

IF(BatchMode) THEN
  StochFile=FirstArg
  nArgsAdd=1
  SWRITE(*,*) "BATCH MODE"
ELSE 
  nArgsAdd=0
  SWRITE(*,*) "SINGLE MODE"
END IF 

IF(nArgsLoc.LT.1+nArgsAdd) CALL Abort(__STAMP__,'Provide Ini File!')
ALLOCATE(ArgsLoc(nArgsLoc-nArgsAdd))
IF(.NOT.BatchMode) ArgsLoc(1)=FirstArg
DO iArg=2,nArgsLoc
  CALL GET_COMMAND_ARGUMENT(iArg,ArgsLoc(iArg-nArgsAdd))
END DO 


IF(BatchMode)THEN
  ! open StochFile, get attributes nGlobalRuns, nParallelRuns
  CALL OpenDataFile(StochFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.,communicatorOpt=MPI_COMM_ACTIVE)
  CALL ReadAttribute(File_ID,'nGlobalRuns',1,IntScalar=nGlobalRuns)
  CALL ReadAttribute(File_ID,'nParallelRuns',1,IntScalar=nParallelRuns)
  CALL CloseDataFile()
ELSE 
  nGlobalRuns=1
  nParallelRuns=1
END IF 

nProcsPerRun = nGlobalProcessors/nParallelRuns
IF(MOD(nGlobalProcessors,nProcsPerRun).NE.0) CALL Abort(__STAMP__,'nProcs has to be a multiple of nProcsPerRun')
iParallelRun = myGlobalRank/nProcsPerRun+1

CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,iParallelRun,myGlobalRank,MPI_COMM_FLEXI,iError) 

!We now allow nGlobalRuns not to fit perfectly.
!IF(MOD(nGlobalRuns,nParallelRuns).NE.0) CALL Abort(__STAMP__,'nGlobalRuns has to be a multiple of nParallelRuns')
nSequentialRuns = (nGlobalRuns-1) / nParallelRuns + 1

! run FLEXI in loop
GlobalStartTime=FLEXITIME()

DO iSequentialRun=1,nSequentialRuns

  iGlobalRun=iParallelRun+nParallelRuns*(iSequentialRun-1)

  ! During last sequential runs, some parallel runs might idle. We therefore split MPI_COMM_WORLD.
  IF(iSequentialRun.EQ.nSequentialRuns)THEN
    isActive=iGlobalRun.LE.nGlobalRuns
    Color=MERGE(0,MPI_UNDEFINED,isActive)
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,Color,myGlobalRank,MPI_COMM_ACTIVE,iError) 
    IF(.NOT.isActive) EXIT
  END IF 

  ! Initialize
  CALL InitFlexi(nArgsLoc-1,ArgsLoc,mpi_comm_loc=MPI_COMM_FLEXI)
  ! Run Simulation
  CALL TimeDisc()

  ! Finalize
  CALL FinalizeFlexi()
END DO
GlobalEndTime=FLEXITIME()

SWRITE(Unit_StdOut,'(A)')
SWRITE(Unit_StdOut,'(A)') "FLEXIBATCH FINISHED"
SWRITE(Unit_StdOut,'(A)') "AvgWork (s) ="
SWRITE(Unit_StdOut,'(E11.5)') (GlobalEndTime-GlobalStartTime)*nGlobalProcessors/(nParallelRuns*nSequentialRuns)

#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif
END PROGRAM FlexiBatch
