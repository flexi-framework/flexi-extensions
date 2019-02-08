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
USE MOD_IO_HDF5    ,ONLY: OpenDataFile,CloseDataFile,File_ID
USE MOD_HDF5_Input ,ONLY: ReadAttribute
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: nArgsLoc,iArg,Color
CHARACTER(LEN=255),ALLOCATABLE :: ArgsLoc(:)
CHARACTER(LEN=255)             :: StochFile
!==================================================================================================================================

! read command line arguments

nArgsLoc = COMMAND_ARGUMENT_COUNT()
IF(nArgsLoc.LT.2) CALL Abort(__STAMP__,'Usage: ./flexi flexi.h5 flexi.ini')
CALL GET_COMMAND_ARGUMENT(1,StochFile)

ALLOCATE(ArgsLoc(nArgsLoc-1))
DO iArg=2,nArgsLoc
  CALL GET_COMMAND_ARGUMENT(iArg,ArgsLoc(iArg-1))
END DO 

! init global MPI

#if USE_MPI
CALL MPI_INIT(iError)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myGlobalRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nGlobalProcessors, iError)

MPI_COMM_ACTIVE=MPI_COMM_WORLD
#else
myGlobalRank=0
nGlobalProcessors=1
#endif

! open StochFile, Get Attributes nRuns, nParallelRuns
CALL OpenDataFile(StochFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL ReadAttribute(File_ID,'nRuns',1,IntScalar=nRuns)
CALL ReadAttribute(File_ID,'nParallelRuns',1,IntScalar=nParallelRuns)
CALL CloseDataFile()


IF(MOD(nGlobalProcessors,nProcsPerRun).NE.0) CALL Abort(__STAMP__,'nProcs has to be a multiple of nProcsPerRun')
nProcsPerRun = nGlobalProcessors/nParallelRuns
iParallelRun = myGlobalRank/nProcsPerRun+1

CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,iParallelRun,myGlobalRank,MPI_COMM_FLEXI,iError) 

!We now allow nRuns not to fit perfectly.
!IF(MOD(nRuns,nParallelRuns).NE.0) CALL Abort(__STAMP__,'nRuns has to be a multiple of nParallelRuns')
nSequentialRuns = nRuns / nParallelRuns

! run FLEXI in loop

DO iSequentialRun=1,nSequentialRuns

  iRun=iSequentialRun+nParallelRuns*(iSequentialRun-1)

  ! During las sequential runs, some parallel runs might idle. We therefore split MPI_COMM_WORLD.
  IF(iSequentialRun.EQ.nSequentialRuns)THEN
    Color=MERGE(0,MPI_UNDEFINED,iRun.LE.nRuns)
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myGlobalRank,MPI_COMM_ACTIVE,iError) 
  END IF 

  ! Initialize
  CALL InitFlexi(nArgsLoc-1,ArgsLoc,mpi_comm_loc=MPI_COMM_FLEXI)
  ! Run Simulation
  CALL TimeDisc()
  ! Finalize
  CALL FinalizeFlexi()
! we also have to finalize MPI itself here
END DO


#if USE_MPI
CALL MPI_FINALIZE(iError)
IF(iError .NE. 0) STOP 'MPI finalize error'
#endif
END PROGRAM FlexiBatch
