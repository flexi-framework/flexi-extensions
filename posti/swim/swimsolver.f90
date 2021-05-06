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
#include "eos.h"

!===================================================================================================================================
!> Containes the routines that will initialize and finalize the SortIcingSides routine as well as the main routine used in the wall
!> distance calculation.
!===================================================================================================================================
MODULE MOD_SwimSolver
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitSwim
  MODULE PROCEDURE InitSwim
END INTERFACE

PUBLIC:: InitSwim

CONTAINS

!===================================================================================================================================
!> Init Swim
!===================================================================================================================================
SUBROUTINE InitSwim()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Swim_Vars
USE MOD_IO_HDF5
USE MOD_HDF5_Input,         ONLY: GetDataSize
USE MOD_HDF5_Input,         ONLY: ReadArray,GetArrayAndName
USE MOD_HDF5_Input,         ONLY: ReadAttribute
USE MOD_HDF5_Output,        ONLY: WriteArray
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iSideIn,SideID,flip,lo,up,iCell,iRun
!===================================================================================================================================


!Read data from state file 
!StateFile   = GETSTR('StateFile')
WRITE (*,*) "Read Data"
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'IceSurfData',nDims,HSize)
CALL ReadAttribute(File_ID,'MeshFile',1,StrScalar=MeshFile)
nVarSurf = INT(HSize(1))
ICS_N = INT(HSize(2))-1
ICS_NZ = INT(HSize(3))-1
nWallSides = INT(HSize(4))
nRuns = INT(HSize(5))
! WriteDim and ICS_NZ should be redundant 
ALLOCATE(VolData(nVarSurf,0:ICS_N,0:ICS_NZ,nWallSides,nRuns))
CALL ReadArray('IceSurfData',5,(/nVarSurf,ICS_N+1,ICS_NZ+1,nWallSides,nRuns/),0,5,RealArray=VolData)
CALL CloseDataFile()


!Read mappingg from mesh file 
!MeshFile   = GETSTR('MeshFile')
WRITE (*,*) "Read Mesh Mapping"
CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataSize(File_ID,'IcingSideInfo',nDims,HSize)
ALLOCATE(Mapping(2,nWallSides))
CALL ReadArray('IcingSideInfo',2,(/2,nWallSides/),0,2,IntArray=Mapping)
CALL CloseDataFile()

nCells = nWallSides*(ICS_N+1)
!Allocate and fill Input arrays
!TODO: in Vars!
ALLOCATE(SwimData(nVarSurf,nCells,nRuns))

WRITE (*,*) "Sort"
DO iSideIn=1,nWallSides
  SideID = Mapping(1,iSideIn)
  lo = (SideID-1)*(ICS_N+1) + 1
  up =  SideID   *(ICS_N+1)
  IF(WriteDim.EQ.3)THEN ! SurfData was created in 3D
    flip   = Mapping(2,iSideIn)
  ELSE ! 2D: Data is already sorted clockwise (due to p SurfVec orientation in 2D CGNS) 
    flip = 1
  END IF 
  Do iRun=1,nRuns
    CALL DoFlip(nVarSurf,VolData(:,:,:,iSideIn,1),flip,SwimData(:,lo:up,iRun))
  END DO 
END DO 

WRITE (*,*) "Write To File"
CALL OpenDataFile(StateFile,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteArray('SwimData',3,(/nVarSurf,nCells,nRuns/),(/nVarSurf,nCells,nRuns/),(/0,0,0/),collective=.FALSE.,RealArray=SwimData)
CALL CloseDataFile()

END SUBROUTINE InitSwim


!===================================================================================================================================
!> Rotate Surf Data
!===================================================================================================================================
SUBROUTINE DoFlip(NIn,UIn,flip,UOut)
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_Swim_Vars
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: NIn
REAL,INTENT(IN)      :: UIn(NIn,0:ICS_N,0:ICS_NZ)
INTEGER,INTENT(IN)   :: flip
REAL,INTENT(OUT)     :: UOut(NIn,0:ICS_N)
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
SELECT CASE(flip)
CASE(1)
  UOut(:,:)=UIn(:,:,0)
CASE(2)
  UOut(:,:)=UIn(:,0,:)
CASE(3)
  UOut(:,:)=UIn(:,ICS_N:0:-1,0)
CASE(4)
  UOut(:,:)=UIn(:,0,ICS_N:0:-1)
END SELECT
END SUBROUTINE DoFlip


END MODULE MOD_SwimSolver
