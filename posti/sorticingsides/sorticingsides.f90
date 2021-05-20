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
MODULE MOD_SortIcingSides
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE CalcSortIcingSides
  MODULE PROCEDURE CalcSortIcingSides
END INTERFACE

PUBLIC:: CalcSortIcingSides

CONTAINS


!===================================================================================================================================
!> Main SortIcingSides routine.
!===================================================================================================================================
SUBROUTINE CalcSortIcingSides()
! MODULES                                                                                                                          !
!----------------------------------------------------------------------------------------------------------------------------------!
USE MOD_PreProc
USE MOD_Globals
USE MOD_ReadInTools
USE MOD_Interpolation_Vars
USE MOD_Mesh_Vars,          ONLY: nBCSides,nGlobalElems,NGeo,NodeCoords
USE MOD_Mesh_Vars,          ONLY: BoundaryType,BC,MeshFile,SideToElem,ElemToSide
USE MOD_Mappings,           ONLY: SideToVol2
USE MOD_HDF5_Output,        ONLY: WriteArray
USE MOD_IO_HDF5
USE MOD_Mesh_ReadIn,        ONLY: ReadIJKSorting
USE MOD_Mesh_Vars,          ONLY: nElems_IJK,Elem_IJK
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER    :: SideID,ElemID,LocSideID,BCType,flip,GlobSideID,nWallSides
INTEGER,ALLOCATABLE    :: SurfIndices(:,:),SurfIndicesAll(:,:)
REAL :: xloc(3,4),xtmp(3,0:1,0:1)

TYPE tSide
  INTEGER                      :: flip
  INTEGER                      :: GlobsideID
  REAL                         :: x(3,4)
  TYPE(tSide),POINTER          :: nextSide
  TYPE(tSide),POINTER          :: prevSide
END TYPE tSide

TYPE tChain
  INTEGER :: i
  TYPE(tSide),POINTER          :: firstSide
  TYPE(tSide),POINTER          :: lastSide
  TYPE(tChain),POINTER         :: nextChain
  TYPE(tChain),POINTER         :: prevChain
END TYPE tChain

TYPE tChainPointer
  TYPE(tChain),POINTER         :: firstChain
END TYPE tChainPointer
TYPE(tChain),POINTER           :: aChain,bChain
TYPE(tChainPointer),ALLOCATABLE:: Chains(:)
TYPE(tSide),POINTER            :: aSide
LOGICAL                        :: found,connected

INTEGER             :: pq(2),p,q
REAL                :: tmp(3,0:NGeo,0:NGeo),xTrailingEdge
INTEGER             :: nZ,nCirc,iZ
!===================================================================================================================================
!clockwise = GETLOGICAL('clockwise','T')
!XCut = GETREALARRAY('XCut',3)

!--------------------------------------------------------------------------------

nElems_IJK = 0 
CALL ReadIJKSorting()
IF(SUM(nElems_IJK).EQ.0) CALL Abort(__STAMP__,'Build Mesh with Elem_IJK!')
nZ = nElems_IJK(3)

ALLOCATE(Chains(nZ))
DO iZ= 1,nZ
  ALLOCATE(Chains(iZ)%firstChain)
  Chains(iZ)%firstChain%i = 1
  NULLIFY(Chains(iZ)%firstChain%firstSide)
  NULLIFY(Chains(iZ)%firstChain%lastSide )
  NULLIFY(Chains(iZ)%firstChain%nextChain)
  NULLIFY(Chains(iZ)%firstChain%prevChain)
END DO 

xTrailingEdge = -1.E10
nCirc = 0 
DO SideID=1,nBCSides
  BCType  = Boundarytype(BC(SideID),BC_TYPE)
  IF(.NOT.ANY(BCType.EQ.WALLBCTYPES())) CYCLE

  ElemID = SideToElem(S2E_ELEM_ID,SideID)
  iZ = Elem_IJK(3,ElemID)

  IF(iZ.EQ.1) nCirc = nCirc + 1 

  ! offsetElem = 0 in single mode!
  locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  GlobSideID = (ElemID-1)*6+LocSideID
  flip = ElemToSide(E2S_FLIP,LocSideID,ElemID)

  SELECT CASE(LocSideID)
  CASE(XI_MINUS)
    tmp=NodeCoords(1:3,0   ,:   ,:   ,ElemID)
  CASE(XI_PLUS)
    tmp=NodeCoords(1:3,NGeo,:   ,:   ,ElemID)
  CASE(ETA_MINUS)
    tmp=NodeCoords(1:3,:   ,0   ,:   ,ElemID)
  CASE(ETA_PLUS)
    tmp=NodeCoords(1:3,:   ,NGeo,:   ,ElemID)
  CASE(ZETA_MINUS)
    tmp=NodeCoords(1:3,:   ,:   ,0   ,ElemID)
  CASE(ZETA_PLUS)
    tmp=NodeCoords(1:3,:   ,:   ,NGeo,ElemID)
  END SELECT
  ! turn into right hand system of side
  DO q=0,ZDIM(NGeo),NGeo; DO p=0,NGeo,NGeo
    pq=SideToVol2(NGeo,p,q,0,LocSideID,PP_dim)
    xtmp(1:3,p/NGeo,q/NGeo)=tmp(:,pq(1),pq(2))
  END DO; END DO ! p,q
  CALL GetFlip(xtmp,xloc,flip)
  xTrailingEdge = MAX(MAXVAL(xloc(1,:)),xTrailingEdge)


  ! TRY TO APPEND TO CHAIN
  aChain => Chains(iZ)%firstChain
  NULLIFY(aSide)
  DO WHILE(ASSOCIATED(aChain))
    found = .FALSE.
    connected = .FALSE.
    IF(ASSOCIATED(aChain%firstSide))THEN
      IF(CONNECT(xloc(:,2),aChain%firstSide%x(:,1)))THEN
        IF(ASSOCIATED(aSide))THEN !Connect bChain --> aChain
          bChain%lastSide%nextSide => aChain%firstSide
          aChain%firstSide%prevSide => bChain%lastSide
          bChain%lastSide => aChain%lastSide
          connected=.TRUE.
        ELSE ! prepend side to chain
          ALLOCATE(aChain%firstSide%prevSide)
          aSide => aChain%firstSide%prevSide
          NULLIFY(aSide%prevSide)
          aSide%nextSide => aChain%firstSide
          aChain%firstSide => aSide
          found = .TRUE.
        ENDIF 
      ELSEIF(CONNECT(xloc(:,1),aChain%lastSide%x(:,2)))THEN
        IF(ASSOCIATED(aSide))THEN !Connect aChain --> bChain
          aChain%lastSide%nextSide => bChain%firstSide
          bChain%firstSide%prevSide => aChain%lastSide
          bChain%firstSide => aChain%firstSide
          connected=.TRUE.
        ELSE ! append side to chain
          ALLOCATE(aChain%lastSide%nextSide)
          aSide => aChain%lastSide%nextSide
          NULLIFY(aSide%nextSide)
          aSide%prevSide => aChain%lastSide
          aChain%lastSide => aSide
          found = .TRUE.
        END IF 
      END IF 
    ELSE 
      IF(.NOT.ASSOCIATED(aSide))THEN
        ! FOUND IN NO PREVIOUS CHAIN: START NEW CHAIN
        ALLOCATE(aChain%firstSide)
        aSide => aChain%firstSide
        NULLIFY(aSide%nextSide)
        NULLIFY(aSide%prevSide)
        aChain%lastSide => aSide
        ALLOCATE(aChain%nextChain)
        aChain%nextChain%i = aChain%i + 1
        NULLIFY(aChain%nextChain%firstSide)
        NULLIFY(aChain%nextChain%lastSide )
        NULLIFY(aChain%nextChain%nextChain)
        !NULLIFY(aChain%nextChain%prevChain)
        aChain%nextChain%prevChain => aChain
        found = .TRUE.
      END IF 
    END IF 
    IF(found)THEN
      aSide%flip = flip
      aSide%x = xloc
      aSide%GlobSideID = GlobSideID
      bChain => aChain
    ELSEIF(connected)THEN
      ! remove aChain and exit
      IF(ASSOCIATED(aChain%nextChain))THEN
        aChain%prevChain%nextChain => aChain%nextChain
        aChain%nextChain%prevChain => aChain%prevChain
      ELSE
        NULLIFY(aChain%prevChain%nextChain)
      ENDIF 
      DEALLOCATE(aChain) 
      EXIT
    END IF 
    aChain => aChain%nextChain
  END DO 
END DO 

! join circle 
DO iZ= 1,nZ
  Chains(iZ)%firstChain%firstSide%prevSide => Chains(iZ)%firstChain%lastSide
  Chains(iZ)%firstChain%lastSide%nextSide => Chains(iZ)%firstChain%firstSide
  !cut circle at back
  aSide => Chains(iZ)%firstChain%firstSide
  DO
    IF(ABS(aSide%x(1,1)-xTrailingEdge).LT.1.E-10) EXIT
    aSide => aSide%nextSide
  END DO 
  Chains(iZ)%firstChain%firstSide => aSide
  NULLIFY(aSide%prevSide%nextSide)
END DO 

!--------------------------------------------------------------------------------
ALLOCATE(SurfIndicesAll(3,6*nGlobalElems))
ALLOCATE(SurfIndices(3,nCirc*nZ))
SurfIndicesAll=0

DO iZ= 1,nZ
  aSide => Chains(iZ)%firstChain%firstSide
  nCirc = 0 
  DO WHILE(ASSOCIATED(aSide))
    nCirc = nCirc + 1
    SurfIndicesAll(1,aSide%GlobSideID) = nCirc
    SurfIndicesAll(2,aSide%GlobSideID) = iZ
    SurfIndicesAll(3,aSide%GlobSideID) = aSide%flip
    aSide => aSide%nextSide
  END DO 
END DO 

! reduce array
nWallSides = 0 
DO GlobSideID = 1,6*nGlobalElems
  IF(SurfIndicesAll(1,GlobSideID).GT.0)THEN
    nWallSides = nWallSides + 1
    SurfIndices(:,nWallSides) = SurfIndicesAll(:,GlobSideID)
  END IF
END DO 

CALL OpenDataFile(MeshFile,create=.FALSE.,single=.FALSE.,readOnly=.FALSE.)
CALL WriteArray("IcingSideInfo",2,(/3,nWallSides/),(/3,nWallSides/),(/0,0/),collective=.TRUE.,IntArray=SurfIndices)
CALL CloseDataFile()

DEALLOCATE(SurfIndicesAll)
DEALLOCATE(SurfIndices)
DEALLOCATE(NodeCoords)

!TODO: Deallocate tSide and tChain objects to free memory
END SUBROUTINE CalcSortIcingSides

!===================================================================================================================================
!> Assumes extrusion in Z.
!> returns p running clockwise, which is true for flip = 1 (which is q in positive Z)
!> this works because q is 90 degree clockwise of p for all iLocSide
!===================================================================================================================================
SUBROUTINE GetFLip(xin,xloc,flip)
! MODULES                                                                                                                          !
USE MOD_Globals
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
REAL,INTENT(IN)      :: xin(3,0:1,0:1)
REAL,INTENT(OUT)     :: xloc(3,4)
INTEGER,INTENT(OUT)  :: flip
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!===================================================================================================================================
xloc(1,:) = (/xin(1,0,0), xin(1,1,0), xin(1,1,1), xin(1,0,1)/)
xloc(2,:) = (/xin(2,0,0), xin(2,1,0), xin(2,1,1), xin(2,0,1)/)
xloc(3,:) = (/xin(3,0,0), xin(3,1,0), xin(3,1,1), xin(3,0,1)/)
flip = 1
DO 
  IF ((xloc(3,3)-xloc(3,2)).GT.1.E-10) RETURN
  !rotate 
  xloc(1,:) = (/xloc(1,2), xloc(1,3), xloc(1,4), xloc(1,1)/)
  xloc(2,:) = (/xloc(2,2), xloc(2,3), xloc(2,4), xloc(2,1)/)
  xloc(3,:) = (/xloc(3,2), xloc(3,3), xloc(3,4), xloc(3,1)/)
  flip = flip + 1
  IF(flip.EQ.5) CALL Abort(__STAMP__,'Error in FLIP')
END DO 
END SUBROUTINE GetFLip

!==================================================================================================================================
!> True if x1 and x2 are equal and not equal to xcut.
!==================================================================================================================================
PPURE FUNCTION CONNECT(x1,x2)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN) :: x1(3)   !< input vector 1
REAL,INTENT(IN) :: x2(3)   !< input vector 2
LOGICAL         :: CONNECT !< cross product of vectors
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CONNECT=((SUM(ABS(x1-x2)).LT.1.E-10))
END FUNCTION CONNECT

END MODULE MOD_SortIcingSides
