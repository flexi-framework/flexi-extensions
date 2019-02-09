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
USE MOD_TimeDisc,          ONLY:TimeDisc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: nParallelRuns=2
INTEGER :: nSequentialRuns,iArg,iSample,nArgsLoc,nProcsPerRun,myParallelRun
CHARACTER(LEN=255),ALLOCATABLE:: ArgsLoc(:) 
!==================================================================================================================================
! Initialize

CALL MPI_INIT(iError)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myGlobalRank     , iError)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nGlobalProcs     , iError)

nGlobalRuns = 4
nArgsLoc = COMMAND_ARGUMENT_COUNT()
!IF(nArgsLoc.LT.2) CALL Abort(__STAMP__,'Usage: ./flexi flexi.h5 flexi.ini')
!CALL GET_COMMAND_ARGUMENT(1,FileNameStochH5)

!! open FileNameStochH5, Get Attributes nGlobalRuns, nParallelRuns

MPI_COMM_ACTIVE=MPI_COMM_WORLD
ALLOCATE(ArgsLoc(nArgsLoc-1))
DO iArg=2,nArgsLoc
  CALL GET_COMMAND_ARGUMENT(iArg,ArgsLoc(iArg-1))
END DO 

!pos_NRuns=MINLOC(


!prmfile
!nGlobalRuns
!nParallelRuns

!IF(MOD(nGlobalProcessors,nProcsPerRun).NE.0) CALL Abort(__STAMP__,'nProcs has to be a multiple of nProcsPerRun')
nProcsPerRun  = nGlobalProcs/nParallelRuns
myParallelRun = myGlobalRank/nProcsPerRun

MPIGlobalRoot=(myGlobalRank .EQ. 0)
CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myParallelRun,myGlobalRank,MPI_COMM_FLEXI,iError) 

!IF(MOD(nRuns,nParallelRuns).NE.0) CALL Abort(__STAMP__,'nRuns has to be a multiple of nParallelRuns')
nSequentialRuns= nGlobalRuns/nParallelRuns
!#else
!nParallelRuns=1
!nSequentialRuns=nGlobalRuns
!#endif

DO iSequentialRun=1,nSequentialRuns
  myGlobalRun = myParallelRun+nParallelRuns*(iSequentialRun-1)
  CALL InitFlexi(nArgsLoc-1,ArgsLoc,MPI_COMM_FLEXI)
  !! Run Simulation
  CALL TimeDisc()
  !! Finalize
  CALL FinalizeFlexi()
#if USE_MPI
  CALL MPI_BARRIER(MPI_COMM_ACTIVE,iError)
#endif
!! we also have to finalize MPI itself here
END DO


!#if USE_MPI
!CALL MPI_FINALIZE(iError)
!IF(iError .NE. 0) STOP 'MPI finalize error'
!#endif
END PROGRAM FlexiBatch
