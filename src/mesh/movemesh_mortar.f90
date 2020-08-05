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

!===================================================================================================================================
!> Contains subroutines necessary for moving meshes in combination with mortar interfaces
!===================================================================================================================================
MODULE MOD_MoveMesh_Mortar
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitMoveMeshMortar
  MODULE PROCEDURE InitMoveMeshMortar
END INTERFACE

INTERFACE MortarMeshmovement
  MODULE PROCEDURE MortarMeshmovement
END INTERFACE

PUBLIC::InitMoveMeshMortar,MortarMeshmovement
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Initialize arrays and communication patterns for mortar interfaces in combination with moving meshes.
!===================================================================================================================================
SUBROUTINE InitMoveMeshMortar()
! MODULES                                                                                                                          !
USE MOD_Globals
USE MOD_PreProc
USE MOD_MoveMesh_Vars
USE MOD_Mesh_Vars,            ONLY: NGeo,MortarType,MortarInfo
USE MOD_Mesh_Vars,            ONLY: firstMPISide_MINE,lastMPISide_MINE,firstMPISide_YOUR,lastMPISide_YOUR
USE MOD_Mesh_Vars,            ONLY: firstInnerSide,firstMortarInnerSide,lastMortarInnerSide,firstMortarMPISide,lastMortarMPISide
USE MOD_Interpolation_Vars,   ONLY: NodeTypeCL
USE MOD_Mortar_Vars,          ONLY: M_0_1_Mesh,M_0_2_Mesh
USE MOD_Mortar,               ONLY: MortarBasis_BigToSmall
#if USE_MPI
USE MOD_MPI_Vars,             ONLY: MPIRequest_MeshMortarRcv,MPIRequest_MeshMortarSend,DataSizeSideMesh,nbProc,nNbProcs
USE MOD_MPI_Vars,             ONLY: OffsetMPISides_rec,OffsetMPISides_send
#endif
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: SideID
INTEGER             :: iNbProc,SendID,iMortar,locSide,MortarSideID,nMortars
INTEGER             :: nSmallMortarMPI_YOURSidesSend,nSmallMortarMPI_MINESidesSend
INTEGER,ALLOCATABLE :: mapSmallMortarsMPI_MINE(:),mapSmallMortarsMPI_YOUR(:)
INTEGER             :: foundNumber,smallestNumber,iSmallMortarsSide,i,j
INTEGER,ALLOCATABLE :: SmallMortarMPISides_tmp(:,:)
INTEGER,ALLOCATABLE :: NbProcsMortar_tmp(:)
INTEGER             :: iSmallMortarMPISides,iSmallMortarSides
INTEGER,ALLOCATABLE :: SmallMortarMPISidesSend(:,:)
!===================================================================================================================================
! Prepare mortar operators for mesh based operations
ALLOCATE(M_0_1_Mesh(0:NGeo,0:NGeo))
ALLOCATE(M_0_2_Mesh(0:NGeo,0:NGeo))
CALL MortarBasis_BigToSmall(0,NGeo,NodeTypeCL, M_0_1_Mesh, M_0_2_Mesh)

! Mortar communication:
! We later need to send the displacement and velocities from the big mortar sides to the small ones
! over the MPI borders. This needs a completly new form of communication since we ONLY need to
! communicate over non-conforming MPI borders, and send from MINE as well as YOUR sides.

!-------------- Prepare mortar recieve -----------------!

! Count the number of small mortar sides:
!   * inner sides where we do not need to communicate
!   * small MPI sides belonging to a small element where we need to recieve
nSmallInnerMortarSides  = 0
nSmallMortarMPISidesRcv = 0
DO SideID=firstInnerSide,lastMPISide_YOUR
  IF (MortarType(1,SideID).EQ.-1) THEN ! This is a small mortar side
    IF (SideID.LT.firstMPISide_MINE) THEN
      nSmallInnerMortarSides  = nSmallInnerMortarSides  + 1
#if USE_MPI
    ELSE
      nSmallMortarMPISidesRcv = nSmallMortarMPISidesRcv + 1
#endif
    END IF
  END IF
END DO

! Build arrays with SideIDs of small mortars, for MPI additionally search for the rank of the corresponding neighbour proc
ALLOCATE(SmallInnerMortarSideIDs( nSmallInnerMortarSides))
#if USE_MPI
ALLOCATE(SmallMortarMPISidesRcv(2,nSmallMortarMPISidesRcv)) ! First Index: SideID, Second Index: iProc
iSmallMortarMPISides = 0
#endif
iSmallMortarSides    = 0
nNbProcsMortarRcv    = 0 ! Count number of procs I recieve mesh movement info's from
DO SideID=firstInnerSide,lastMPISide_YOUR
  IF (MortarType(1,SideID).EQ.-1) THEN ! This is a small mortar side
    IF (SideID.LT.firstMPISide_MINE) THEN
      ! inner side
      iSmallMortarSides = iSmallMortarSides + 1
      SmallInnerMortarSideIDs(iSmallMortarSides) = SideID
#if USE_MPI
    ELSE IF (SideID.LE.lastMPISide_YOUR) THEN
      ! MPI MINE or YOUR side
      IF (SideID.GE.firstMPISide_YOUR) THEN
        SendID=2
      ELSE
        SendID=1
      END IF
      iSmallMortarMPISides = iSmallMortarMPISides + 1
      SmallMortarMPISidesRcv(1,iSmallMortarMPISides) = SideID
      DO iNbProc = 1, nNbProcs
        IF ((SideID.GT.OffsetMPISides_send(iNbProc-1,SendID)).AND.(SideID.LE.OffsetMPISides_send(iNbProc,SendID))) THEN
          SmallMortarMPISidesRcv(2,iSmallMortarMPISides) = nbProc(iNbProc)
          EXIT
        END IF
      END DO ! iNbProc = 1, nNbProcs
      ! Increase count of rcv neighbour procs if this rank has not yet been added
      IF (.NOT.ANY(SmallMortarMPISidesRcv(2,1:iSmallMortarMPISides-1).EQ.SmallMortarMPISidesRcv(2,iSmallMortarMPISides))) THEN
        nNbProcsMortarRcv = nNbProcsMortarRcv + 1
      END IF
#endif
    END IF
  END IF
END DO

#if USE_MPI
! Build array of ranks I recieve from
ALLOCATE(NbProcsMortarRcv(nNbProcsMortarRcv))
iNbProc = 0
DO SideID=1,nSmallMortarMPISidesRcv
  IF (.NOT.ANY(SmallMortarMPISidesRcv(2,1:SideID-1).EQ.SmallMortarMPISidesRcv(2,SideID))) THEN
    ! New rank
    iNbProc = iNbProc + 1
    NbProcsMortarRcv(iNbProc) = SmallMortarMPISidesRcv(2,SideID)
  END IF
END DO
! Order them in ascending order
ALLOCATE(NbProcsMortar_tmp(nNbProcsMortarRcv))
foundNumber = -1
iNbProc = 0
DO i = 1, nNbProcsMortarRcv
  smallestNumber = HUGE(1)
  DO j = 1, nNbProcsMortarRcv
    IF ((NbProcsMortarRcv(j).LT.smallestNumber).AND.((NbProcsMortarRcv(j).GT.foundNumber))) THEN
      smallestNumber  = NbProcsMortarRcv(j)
    END IF
  END DO ! j = 1, nNbProcsMortarRcv
  iNbProc = iNbProc + 1
  NbProcsMortar_tmp(iNbProc) = smallestNumber
  foundNumber = smallestNumber
END DO ! i = 1, nNbProcsMortarRcv
NbProcsMortarRcv = NbProcsMortar_tmp
DEALLOCATE(NbProcsMortar_tmp)

! Count number of sides per proc
ALLOCATE(nMortarSidesRcv(nNbProcsMortarRcv))
nMortarSidesRcv = 0
DO iNbProc=1,nNbProcsMortarRcv
  nMortarSidesRcv(iNbProc) = COUNT(SmallMortarMPISidesRcv(2,:).EQ.NbProcsMortarRcv(iNbProc))
END DO
! Build offset array
ALLOCATE(offsetMortarSidesRcv(nNbProcsMortarRcv+1))
offsetMortarSidesRcv(1) = 0
DO iNbProc=1,nNbProcsMortarRcv
  offsetMortarSidesRcv(iNbProc+1) = offsetMortarSidesRcv(iNbProc) + nMortarSidesRcv(iNbProc)
END DO
! Build map from recieve array to real side IDs - sorted first after the rank and then after the SideID
ALLOCATE(mapMortarSidesRcv(nSmallMortarMPISidesRcv))
i = 0
DO iNbProc=1,nNbProcsMortarRcv
  DO SideID=1,nSmallMortarMPISidesRcv
    IF (SmallMortarMPISidesRcv(2,SideID).EQ.NbProcsMortarRcv(iNbProc)) THEN
      i = i+1
      mapMortarSidesRcv(i) = SmallMortarMPISidesRcv(1,SideID)
    END IF
  END DO
END DO
#endif

!----------------- Prepare mortar send --------------------!
! Loop over all big mortars to find small sides which need to be send
! We don't just loop over all MPI mortars since big mortars with only MPI_MINE small sides are 
! included in the mortarInner sides!
nSmallMortarMPISidesSend = 0
#if USE_MPI
DO MortarSideID=firstMortarInnerSide,lastMortarMPISide
  IF ((MortarSideID.GT.lastMortarInnerSide).AND.(MortarSideID.LT.firstMortarMPISide)) CYCLE
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID) ! ID in mortar info array (1..nBigMortarSides)
  DO iMortar=1,nMortars ! Loop over all small mortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide) ! SideID of small mortar side
    IF (SideID.GE.firstMPISide_MINE) THEN ! This small sides needs to be communicated
      nSmallMortarMPISidesSend = nSmallMortarMPISidesSend + 1
    END IF
  END DO
END DO

! Allocate arrays and fill with send information
ALLOCATE(SmallMortarMPISidesSend(2,nSmallMortarMPISidesSend)) ! First Index: SideID, Second Index: iProc
iSmallMortarMPISides = 0
nNbProcsMortarSend = 0
nSmallMortarMPI_MINESidesSend = 0
nSmallMortarMPI_YOURSidesSend = 0
DO MortarSideID=firstMortarInnerSide,lastMortarMPISide
  IF ((MortarSideID.GT.lastMortarInnerSide).AND.(MortarSideID.LT.firstMortarMPISide)) CYCLE
  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID) ! ID in mortar info array (1..nBigMortarSides)
  DO iMortar=1,nMortars ! Loop over all small mortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide) ! SideID of small mortar side
    IF (SideID.GE.firstMPISide_MINE) THEN ! This small sides needs to be communicated
      iSmallMortarMPISides = iSmallMortarMPISides + 1
      SmallMortarMPISidesSend(1,iSmallMortarMPISides) = SideID
      ! MPI MINE or YOUR side
      IF (SideID.GE.firstMPISide_YOUR) THEN
        nSmallMortarMPI_MINESidesSend = nSmallMortarMPI_MINESidesSend + 1
        SendID=1
      ELSE
        nSmallMortarMPI_YOURSidesSend = nSmallMortarMPI_YOURSidesSend + 1
        SendID=2
      END IF
      DO iNbProc = 1, nNbProcs
        IF ((SideID.GT.OffsetMPISides_rec(iNbProc-1,SendID)).AND.(SideID.LE.OffsetMPISides_rec(iNbProc,SendID))) THEN
          SmallMortarMPISidesSend(2,iSmallMortarMPISides) = nbProc(iNbProc)
          EXIT
        END IF
      END DO ! iNbProc = 1, nNbProcs
      IF (.NOT.(ANY(SmallMortarMPISidesSend(2,1:iSmallMortarMPISides-1).EQ.SmallMortarMPISidesSend(2,iSmallMortarMPISides)))) THEN
        nNbProcsMortarSend = nNbProcsMortarSend + 1
      END IF
    END IF
  END DO
END DO

! Build array of ranks I send to
ALLOCATE(NbProcsMortarSend(nNbProcsMortarSend))
iNbProc = 0
DO SideID=1,nSmallMortarMPISidesSend
  IF (.NOT.(ANY(SmallMortarMPISidesSend(2,1:SideID-1).EQ.SmallMortarMPISidesSend(2,SideID)))) THEN
    ! New rank
    iNbProc = iNbProc + 1
    NbProcsMortarSend(iNbProc) = SmallMortarMPISidesSend(2,SideID)
  END IF
END DO
! Order them in ascending order
ALLOCATE(NbProcsMortar_tmp(nNbProcsMortarSend))
foundNumber = -1
iNbProc = 0
DO i = 1, nNbProcsMortarSend
  smallestNumber = HUGE(1)
  DO j = 1, nNbProcsMortarSend
    IF ((NbProcsMortarSend(j).LT.smallestNumber).AND.((NbProcsMortarSend(j).GT.foundNumber))) THEN
      smallestNumber  = NbProcsMortarSend(j)
    END IF
  END DO ! j = 1, nNbProcsMortarSend
  iNbProc = iNbProc + 1
  NbProcsMortar_tmp(iNbProc) = smallestNumber
  foundNumber = smallestNumber
END DO ! i = 1, nNbProcsMortarSend
NbProcsMortarSend = NbProcsMortar_tmp
DEALLOCATE(NbProcsMortar_tmp)

! The array SmallMortarMPISidesSend contains the sides in more or less "random" fashion. 
! We need to sort the sides and also switch around the MPI_MINE and MPI_YOUR sides to get consistent ordering 
! with the recieving procs (their MINE is my YOUR).
ALLOCATE(mapSmallMortarsMPI_YOUR(nSmallMortarMPI_MINESidesSend))
ALLOCATE(mapSmallMortarsMPI_MINE(nSmallMortarMPI_YOURSidesSend))
foundNumber = -1
DO i = 1, nSmallMortarMPI_YOURSidesSend
  smallestNumber = HUGE(1)
  DO j = 1, nSmallMortarMPISidesSend
    IF (SmallMortarMPISidesSend(1,j).GT.lastMPISide_MINE) CYCLE ! Only for MINE sides
    IF ((SmallMortarMPISidesSend(1,j).LT.smallestNumber).AND.((SmallMortarMPISidesSend(1,j).GT.foundNumber))) THEN
      smallestNumber = SmallMortarMPISidesSend(1,j)
      iSmallMortarsSide = j
    END IF
  END DO ! j = 1, nSmallMortarMPISidesSend
  mapSmallMortarsMPI_MINE(i) = iSmallMortarsSide
  foundNumber = smallestNumber
END DO ! i = 1, nSmallMortarMPI_YOURSidesSend
foundNumber = -1
DO i = 1, nSmallMortarMPI_MINESidesSend
  smallestNumber = HUGE(1)
  DO j = 1, nSmallMortarMPISidesSend
    IF (SmallMortarMPISidesSend(1,j).LT.firstMPISide_YOUR) CYCLE ! Only for YOUR sides
    IF ((SmallMortarMPISidesSend(1,j).LT.smallestNumber).AND.((SmallMortarMPISidesSend(1,j).GT.foundNumber))) THEN
      smallestNumber  = SmallMortarMPISidesSend(1,j)
      iSmallMortarsSide = j
    END IF
  END DO ! j = 1, nSmallMortarMPISidesSend
  mapSmallMortarsMPI_YOUR(i) = iSmallMortarsSide
  foundNumber = smallestNumber
END DO ! i = 1, nSmallMortarMPI_MINESidesSend

! Resort the array SmallMortarMPISidesSend
ALLOCATE(SmallMortarMPISides_tmp(2,nSmallMortarMPISidesSend))
DO iSmallMortarsSide = 1,nSmallMortarMPI_MINESidesSend
  SmallMortarMPISides_tmp(:,iSmallMortarsSide) = SmallMortarMPISidesSend(:,mapSmallMortarsMPI_YOUR(iSmallMortarsSide))
END DO
DO iSmallMortarsSide = 1,nSmallMortarMPI_YOURSidesSend
  SmallMortarMPISides_tmp(:,iSmallMortarsSide+nSmallMortarMPI_MINESidesSend) = &
                         SmallMortarMPISidesSend(:,mapSmallMortarsMPI_MINE(iSmallMortarsSide))
END DO
SmallMortarMPISidesSend = SmallMortarMPISides_tmp
DEALLOCATE(SmallMortarMPISides_tmp)

! Count number of sides per proc
ALLOCATE(nMortarSidesSend(nNbProcsMortarSend))
nMortarSidesSend = 0
DO iNbProc=1,nNbProcsMortarSend
  nMortarSidesSend(iNbProc) = COUNT(SmallMortarMPISidesSend(2,:).EQ.NbProcsMortarSend(iNbProc))
END DO
! Build offset array
ALLOCATE(offsetMortarSidesSend(nNbProcsMortarSend+1))
offsetMortarSidesSend(1) = 0
DO iNbProc=1,nNbProcsMortarSend
  offsetMortarSidesSend(iNbProc+1) = offsetMortarSidesSend(iNbProc) + nMortarSidesSend(iNbProc)
END DO
! Build map from recieve array to real side IDs - sorted first after the rank and then after the SideID
ALLOCATE(mapMortarSidesSend(nSmallMortarMPISidesSend))
i = 0
DO iNbProc=1,nNbProcsMortarSend
  DO SideID=1,nSmallMortarMPISidesSend
    IF (SmallMortarMPISidesSend(2,SideID).EQ.NbProcsMortarSend(iNbProc)) THEN
      i = i+1
      mapMortarSidesSend(i) = SmallMortarMPISidesSend(1,SideID)
    END IF
  END DO
END DO

! Necessary MPI variables
DataSizeSideMesh = 3*(NGeo+1)*(ZDIM(NGeo)+1)
ALLOCATE(MPIRequest_MeshMortarSend(2,nNbProcsMortarSend))
ALLOCATE(MPIRequest_MeshMortarRcv( 2,nNbProcsMortarRcv))
MPIRequest_MeshMortarRcv  = MPI_REQUEST_NULL
MPIRequest_MeshMortarSend = MPI_REQUEST_NULL

! Cleanup
DEALLOCATE(SmallMortarMPISidesSend)
DEALLOCATE(mapSmallMortarsMPI_MINE)
DEALLOCATE(mapSmallMortarsMPI_YOUR)
#endif

END SUBROUTINE InitMoveMeshMortar

!===================================================================================================================================
!> For mortar meshes, we need to make sure that we retain a continuous mesh representation over the mortar interfaces. This will not
!> be automatically the case if the mesh is not able to exactly represent the movement. To ensure this, we need to interpolate the
!> mesh movement from the big mortar sides to the small mortar sides.
!===================================================================================================================================
SUBROUTINE MortarMeshmovement(Displacement,Velocity)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: NGeo,nElems,nSides
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_FillMortarMesh   ,ONLY: Mesh_Mortar
USE MOD_MoveMesh_Vars    ,ONLY: v_NGeo_Face
#if USE_MPI
USE MOD_MoveMesh_Vars    ,ONLY: nSmallMortarMPISidesSend,mapMortarSidesRcv,mapMortarSidesSend
USE MOD_MoveMesh_Vars    ,ONLY: nSmallMortarMPISidesRcv,nNbProcsMortarRcv,nNbProcsMortarSend
USE MOD_MPI_Vars         ,ONLY: MPIRequest_MeshMortarRcv,MPIRequest_MeshMortarSend
USE MOD_MPI              ,ONLY: FinishExchangeMPIData,StartExchange_MeshMortar
#endif
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)    :: Displacement(3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh displacement in volume
REAL,INTENT(INOUT)    :: Velocity(    3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh velocity     in volume
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL               :: Displacement_Face(1:3,0:NGeo,0:ZDIM(NGeo),1:nSides)
#if USE_MPI
REAL               :: SendBufferDisp(1:3,0:NGeo,0:ZDIM(NGeo),1:nSmallMortarMPISidesSend)
REAL               :: RcvBufferDisp( 1:3,0:NGeo,0:ZDIM(NGeo),1:nSmallMortarMPISidesRcv)
REAL               :: SendBufferVel( 1:3,0:NGeo,0:ZDIM(NGeo),1:nSmallMortarMPISidesSend)
REAL               :: RcvBufferVel(  1:3,0:NGeo,0:ZDIM(NGeo),1:nSmallMortarMPISidesRcv)
INTEGER            :: iSide
#endif
!===================================================================================================================================
! General structure:
! 1.): Get the mesh movement on the big mortar sides from the volume data (since we use CL points, the surface data
! directly comes from the volume data by taking a slice)
! 2.): "Mortarize" the quantities -> interpolate them to the small sides (in the right hand system of the small mortar
! side)
! 3.) Communicate the mesh movement from big Mortars to small mortars if necessary (MPI only)
! 4.) Sort the exchanged surface information back into the volume where it is later needed

#if USE_MPI
! 1.) and 2.) for sides that need to be communicated
CALL GetMortarFaceData(Displacement,Velocity,Displacement_Face,v_NGeo_Face,doMPISides=.TRUE.)
CALL Mesh_Mortar(Displacement_Face,doMPIsides=.TRUE.)
CALL Mesh_Mortar(v_NGeo_Face,      doMPIsides=.TRUE.)
#endif
! 1.) and 2.) for sides that are local to this proc
CALL GetMortarFaceData(Displacement,Velocity,Displacement_Face,v_NGeo_Face,doMPISides=.FALSE.)
CALL Mesh_Mortar(Displacement_Face,doMPIsides=.FALSE.)
CALL Mesh_Mortar(v_NGeo_Face,      doMPIsides=.FALSE.)

#if USE_MPI
! 3.)
! Fill the send buffer
DO iSide = 1, nSmallMortarMPISidesSend
  SendBufferDisp(:,:,:,iSide) = Displacement_Face(:,:,:,mapMortarSidesSend(iSide))
  SendBufferVel( :,:,:,iSide) = v_NGeo_Face(      :,:,:,mapMortarSidesSend(iSide))
END DO ! iSide = 1, nSmallMortarMPISidesSend

CALL StartExchange_MeshMortar(SendBufferDisp,RcvBufferDisp,SendBufferVel,RcvBufferVel)
#endif

!4.)
CALL MortarFaceToVolume(Displacement,Velocity,Displacement_Face,v_NGeo_Face,doMPISides=.FALSE.)

#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcsMortarSend,MPIRequest_MeshMortarSend)
CALL FinishExchangeMPIData(2*nNbProcsMortarRcv, MPIRequest_MeshMortarRcv)

! Sort the recieve buffer back into the global array
DO iSide = 1, nSmallMortarMPISidesRcv
  Displacement_Face(:,:,:,mapMortarSidesRcv(iSide)) = RcvBufferDisp(:,:,:,iSide)
  v_NGeo_Face(      :,:,:,mapMortarSidesRcv(iSide)) = RcvBufferVel( :,:,:,iSide)
END DO ! iSide = 1, nSmallMortarMPISidesRcv

CALL MortarFaceToVolume(Displacement,Velocity,Displacement_Face,v_NGeo_Face,doMPISides=.TRUE.)
#endif

END SUBROUTINE MortarMeshmovement

!===================================================================================================================================
!> Use the face data that has been interpolated from the big mortars to the small ones to fill the volume information
!> that is later needed. After this step, we have a continuous mesh representation over the mortar faces.
!===================================================================================================================================
SUBROUTINE MortarFaceToVolume(Displacement,Velocity,Displacement_Face,Velocity_Face,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: ElemToSide,NGeo,SideToElem,nElems,nSides
#if FV_ENABLED
USE MOD_MoveMesh_Vars    ,ONLY: Vdm_CLNGeo_FV
USE MOD_FV_Vars          ,ONLY: FV_Elems_sum
#endif
USE MOD_ChangeBasisByDim ,ONLY: ChangeBasisSurf
USE MOD_FillMortarMesh   ,ONLY: Mesh_Mortar
USE MOD_MoveMesh_Vars    ,ONLY: S2V2_NGeo,nSmallInnerMortarSides,SmallInnerMortarSideIDs
#if USE_MPI
USE MOD_MoveMesh_Vars    ,ONLY: nSmallMortarMPISidesRcv,SmallMortarMPISidesRcv
#endif
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)    :: Displacement(3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh displacement in volume
REAL,INTENT(INOUT)    :: Velocity(    3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh velocity     in volume
REAL,INTENT(INOUT)    :: Displacement_Face(1:3,0:NGeo,0:ZDIM(NGeo),1:nSides)!< Mesh displacement on face
REAL,INTENT(INOUT)    :: Velocity_Face(    1:3,0:NGeo,0:ZDIM(NGeo),1:nSides)!< Mesh velocity     on face
LOGICAL,INTENT(IN)    :: doMPISides                                       !< Switch to do only or no MPI sides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2),iSide
INTEGER            :: ElemID,locSideID,SideID
INTEGER            :: flip,nSidesLoc
INTEGER,POINTER    :: SideIDs(:)
!===================================================================================================================================
#if USE_MPI
IF (doMPISides) THEN
  nSidesLoc =  nSmallMortarMPISidesRcv
  SideIDs   => SmallMortarMPISidesRcv(1,:)
ELSE
#endif
  nSidesLoc =  nSmallInnerMortarSides
  SideIDs   => SmallInnerMortarSideIDs
#if USE_MPI
END IF
#endif

! Loop over all small sides belonging to a big mortar
DO iSide=1,nSidesLoc
  SideID= SideIDs(iSide)
  ElemID= SideToElem(S2E_ELEM_ID,SideID)
  IF (ElemID.LE.0) THEN
    ! The small element is slave to the small side
    ElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)
    locSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
  ELSE
    ! The small element is master to the small side
    locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
  END IF
  ! Flip from small element side to master system on the small mortar side (the system the mesh quantites are in)
  flip = ElemToSide(E2S_FLIP,locSideID,ElemID)

  ! Sort back into volume
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,0    ,pq(1),pq(2),ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,0    ,pq(1),pq(2),ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  CASE(XI_PLUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,NGeo ,pq(1),pq(2),ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,NGeo ,pq(1),pq(2),ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  CASE(ETA_MINUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,pq(1),0    ,pq(2),ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,pq(1),0    ,pq(2),ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  CASE(ETA_PLUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,pq(1),NGeo ,pq(2),ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,pq(1),NGeo ,pq(2),ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  CASE(ZETA_MINUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,pq(1),pq(2),0    ,ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,pq(1),pq(2),0    ,ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  CASE(ZETA_PLUS)
    DO q=0,ZDIM(NGeo); DO p=0,NGeo
      pq=S2V2_NGeo(:,p,q,flip,locSideID)
      Displacement(1:3,pq(1),pq(2),NGeo ,ElemID) = Displacement_Face(1:3,p,q,SideID)
      Velocity(    1:3,pq(1),pq(2),NGeo ,ElemID) = Velocity_Face(    1:3,p,q,SideID)
    END DO; END DO ! p,q
  END SELECT
END DO !iSide

END SUBROUTINE MortarFaceToVolume

!===================================================================================================================================
!> Fill face mesh displacement and velocites for big mortar sides by taking a slice from the volume data.
!===================================================================================================================================
SUBROUTINE GetMortarFaceData(Displacement,Velocity,Displacement_Face,Velocity_Face,doMPISides)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars        ,ONLY: ElemToSide,NGeo,SideToElem,nElems,nSides
USE MOD_Mesh_Vars        ,ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars        ,ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_MoveMesh_Vars    ,ONLY: S2V2_NGeo
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(INOUT)    :: Displacement(3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh displacement in volume
REAL,INTENT(INOUT)    :: Velocity(    3,0:NGeo,0:NGeo,0:ZDIM(NGeo),nElems)  !< Mesh velocity     in volume
REAL,INTENT(INOUT)    :: Displacement_Face(1:3,0:NGeo,0:ZDIM(NGeo),1:nSides)!< Mesh displacement on face
REAL,INTENT(INOUT)    :: Velocity_Face(    1:3,0:NGeo,0:ZDIM(NGeo),1:nSides)!< Mesh velocity     on face
LOGICAL,INTENT(IN)    :: doMPISides                                       !< Switch to do only or no MPI sides
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: p,q,pq(2)
INTEGER            :: ElemID,locSideID,MortarSideID
INTEGER            :: firstMortarSideID,lastMortarSideID
INTEGER            :: flip
REAL               :: tmp( 3,0:NGeo,0:ZDIM(NGeo))
REAL               :: tmp2(3,0:NGeo,0:ZDIM(NGeo))
!===================================================================================================================================
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

! Loop over all big mortar sides
DO MortarSideID = firstMortarSideID, lastMortarSideID
  ! The big sides are always master, we need to get the mesh movement from them
  ElemID    = SideToElem(S2E_ELEM_ID,    MortarSideID)
  locSideID = SideToElem(S2E_LOC_SIDE_ID,MortarSideID)
  flip      = ElemToSide(E2S_FLIP,locSideID,ElemID)
  ! Select slice of volume data depending on the local side
  SELECT CASE(locSideID)
  CASE(XI_MINUS)
    tmp =Displacement(1:3,0   ,:   ,:   ,ElemID)
    tmp2=Velocity(    1:3,0   ,:   ,:   ,ElemID)
  CASE(XI_PLUS)
    tmp =Displacement(1:3,NGeo,:   ,:   ,ElemID)
    tmp2=Velocity(    1:3,NGeo,:   ,:   ,ElemID)
  CASE(ETA_MINUS)
    tmp =Displacement(1:3,:   ,0   ,:   ,ElemID)
    tmp2=Velocity(    1:3,:   ,0   ,:   ,ElemID)
  CASE(ETA_PLUS)
    tmp =Displacement(1:3,:   ,NGeo,:   ,ElemID)
    tmp2=Velocity(    1:3,:   ,NGeo,:   ,ElemID)
  CASE(ZETA_MINUS)
    tmp =Displacement(1:3,:   ,:   ,0   ,ElemID)
    tmp2=Velocity(    1:3,:   ,:   ,0   ,ElemID)
  CASE(ZETA_PLUS)
    tmp =Displacement(1:3,:   ,:   ,NGeo,ElemID)
    tmp2=Velocity(    1:3,:   ,:   ,NGeo,ElemID)
  END SELECT

  ! turn into right hand system of big mortar side
  DO q=0,ZDIM(NGeo); DO p=0,NGeo
    pq=S2V2_NGeo(:,p,q,flip,locSideID)
    Displacement_Face(1:3,p,q,MortarSideID)=tmp( :,pq(1),pq(2))
    Velocity_Face(    1:3,p,q,MortarSideID)=tmp2(:,pq(1),pq(2))
  END DO; END DO ! p,q
END DO ! MortarSideID = firstMortarSideID, lastInnerMortarSideID

END SUBROUTINE GetMortarFaceData

END MODULE MOD_MoveMesh_Mortar
