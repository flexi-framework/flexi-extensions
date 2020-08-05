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
!>
!>
!==================================================================================================================================
MODULE MOD_SM
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER,PARAMETER :: S2SM_IPROC      = 1
INTEGER,PARAMETER :: S2SM_IINTERFACE = 2
INTEGER,PARAMETER :: S2SM_IAZ        = 3
INTEGER,PARAMETER :: S2SM_IL         = 4
INTEGER,PARAMETER :: S2SM_IMORTAR    = 5
INTEGER,PARAMETER :: S2SM_ISIDE      = 6
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE DefineParametersSM
  MODULE PROCEDURE DefineParametersSM
END INTERFACE

INTERFACE InitSM
  MODULE PROCEDURE InitSM
END INTERFACE

INTERFACE PrepareSM
  MODULE PROCEDURE PrepareSM
END INTERFACE

INTERFACE GetSMFlip
  MODULE PROCEDURE GetSMFlip
END INTERFACE

INTERFACE FinalizeSM
  MODULE PROCEDURE FinalizeSM
END INTERFACE

PUBLIC::DefineParametersSM,InitSM,PrepareSM,FinalizeSM,GetSMFlip

!==================================================================================================================================

CONTAINS


!==================================================================================================================================
!> Define parameters used by sliding mesh routines
!==================================================================================================================================
SUBROUTINE DefineParametersSM()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
CALL prms%SetSection("SM")
!CALL prms%CreateRealOption('WeightRot', "element weighting of rotor elements for load balancing", '1.')
!==================================================================================================================================
END SUBROUTINE DefineParametersSM



!==================================================================================================================================
!> Basic SM initialization.
!==================================================================================================================================
SUBROUTINE InitSM()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_SM_Vars
USE MOD_Mesh_Vars, ONLY:firstSMSide,lastSMSide,nSMSides,IAmAStatProc,IAmARotProc
#if FV_ENABLED
USE MOD_Mesh_Vars, ONLY:nElems,SideToElem
#endif
IMPLICIT NONE
INTEGER :: iSide
#if FV_ENABLED
INTEGER :: iElem
#endif
!----------------------------------------------------------------------------------------------------------------------------------
IF(SMInitIsDone) CALL CollectiveStop(__STAMP__,'InitSM already called.')

#if FV_ENABLED
ALLOCATE(SM_Elems(1:nElems))
SM_Elems = 0
! Save which Elems have SM Sides (needed to force them to stay DG)
DO iSide=firstSMSide,lastSMSide
  iElem  = SideToElem(S2E_ELEM_ID   ,iSide)
  IF (iElem.GT.0) SM_Elems(iElem) = 1
  iElem  = SideToElem(S2E_NB_ELEM_ID,iSide)
  IF (iElem.GT.0) SM_Elems(iElem) = 1
END DO
#endif /* FV_ENABLED */

ALLOCATE(M_SM_0_1(0:PP_N,0:PP_N,1:nSlidingMeshInterfaces))
ALLOCATE(M_SM_0_2(0:PP_N,0:PP_N,1:nSlidingMeshInterfaces))
ALLOCATE(M_SM_1_0(0:PP_N,0:PP_N,1:nSlidingMeshInterfaces))
ALLOCATE(M_SM_2_0(0:PP_N,0:PP_N,1:nSlidingMeshInterfaces))

ALLOCATE(NormVec_SM (1:3,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(TangVec1_SM(1:3,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(TangVec2_SM(1:3,0:PP_N,0:PP_NZ,1:2*nSMSides))
NormVec_SM  = 0.
TangVec1_SM = 0.
TangVec2_SM = 0.

ALLOCATE(   U_MorStat( 1:PP_nVar,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(UPrim_MorStat(1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(Flux_MorStat( 1:PP_nVar,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(    U_MorRot( 1:PP_nVar,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(UPrim_MorRot( 1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE( Flux_MorRot( 1:PP_nVar,0:PP_N,0:PP_NZ,1:2*nSMSides))

#if PARABOLIC
ALLOCATE(gradUx_MorStat(1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(gradUy_MorStat(1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(gradUz_MorStat(1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(gradUx_MorRot (1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(gradUy_MorRot (1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE(gradUz_MorRot (1:PP_nVarPrim,0:PP_N,0:PP_NZ,1:2*nSMSides))
#if EDDYVISCOSITY
ALLOCATE( muSGS_MorStat(1            ,0:PP_N,0:PP_NZ,1:2*nSMSides))
ALLOCATE( muSGS_MorRot (1            ,0:PP_N,0:PP_NZ,1:2*nSMSides))
#endif /*EDDYVISCOSITY*/
#endif /*PARABOLIC*/

ALLOCATE(SMFlip(firstSMSide:lastSMSide))
DO iSide=firstSMSide,lastSMSide
  CALL GetSMFlip(iSide,SMFlip(iSide))
END DO

ALLOCATE(SM_nMPIMortars_Proc(0:nProcessors-1))
SM_nMPIMortars_Proc = 0
IF(IAmAStatProc) ALLOCATE(SMMortarToStatSide(2,firstSMSide:lastSMSide))
IF(IAmARotProc)  ALLOCATE(SMMortarToRotSide (2,firstSMSide:lastSMSide))
ALLOCATE(SMn_prev(1:nSlidingMeshInterfaces))
SMn_prev=-1

ALLOCATE(Face_vGPSM(1:3,0:PP_N,0:PP_NZ))
Face_vGPSM=0.


SMInitIsDone=.TRUE.
!==================================================================================================================================
END SUBROUTINE InitSM



!==================================================================================================================================
!> Updates mortar operators and neighbour connectivity for the sliding mesh interface.
!> Interpolation and projections matrices are updated in every RK stage (according to analytical rotation angle).
!> SM mortars are sorted (in this hierarchy) by neighbour proc, azimuth, layer and iMortar (1/2).
!> An index matrix from (proc local) mortar sorting to proc local side ID is written. This matrix is only updated when neighbour
!> connections flip, i.e. after each rotor side has passed a full stator side.
!==================================================================================================================================
SUBROUTINE PrepareSM(t)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_SM_Vars
USE MOD_Interpolation_Vars, ONLY:NodeType
USE MOD_Mortar,             ONLY:MortarBasis_SmallToBig,MortarBasis_BigToSmall
USE MOD_MortarSM,           ONLY:DoFlip3,InterpolateSM3
USE MOD_Mesh_Vars,          ONLY:firstSMSide,lastSMSide,SlidingMeshInfo
USE MOD_Mesh_Vars,          ONLY:nAzimuthalSides,nSMSides
USE MOD_Mesh_Vars,          ONLY:IAmAStatProc,IAmARotProc,AL_ToStatProc,AL_ToRotProc
USE MOD_Mesh_Vars,          ONLY:LocToGlobInterface,IntToPart
USE MOD_Mesh_Vars,          ONLY:NormVec,TangVec1,TangVec2
USE MOD_MoveMesh_Vars,      ONLY:RotationAngVel,MeshVel,MeshAcc,DoMeshAccelerate
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(IN)                 :: t                      !< Current time
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,POINTER      :: M1(:,:),M2(:,:)
REAL              :: delta_phi,tmp
REAL              :: VecTmp(1:3,0:PP_N,0:PP_N)
REAL              :: SMo(nSlidingMeshInterfaces)
INTEGER           :: SMn(nSlidingMeshInterfaces)
INTEGER           :: iSide,iPartition,i
INTEGER           :: iAz,iL,iMortar,iProc,iiMortar,iInterface,iGlobInterface
INTEGER           :: Mortars(2)
INTEGER           :: SideToSMMortar(2*nSMSides,6)
INTEGER,PARAMETER :: STAT=1,ROT=2
!==================================================================================================================================

DO iInterface=1,nSlidingMeshInterfaces
  iGlobInterface = LocToGlobInterface(iInterface)
  iPartition     = IntToPart(iGlobInterface)

  ! SOME ANALYTICAL VALUES CALCULATED FROM ROTATION ANGLE
  !----------------------------------------------------------------------------------------------------
  SELECT CASE (SlidingMeshType(iPartition))
    CASE(SM_TYPE_ROTATION_RADIAL,SM_TYPE_ROTATION_AXIAL)
      delta_phi=2.*PP_Pi/REAL(nAzimuthalSides(iGlobInterface))
      tmp=t*RotationAngVel(iPartition) / delta_phi
    CASE(SM_TYPE_PLANAR)
      delta_phi=(SlidingMeshBoundaries(2,iPartition)-SlidingMeshBoundaries(1,iPartition))/REAL(nAzimuthalSides(iGlobInterface))
      IF (DoMeshAccelerate(iPartition)) THEN
        tmp=0.5*t**2*MeshAcc(SlidingMeshDirection(iPartition),iPartition)/delta_phi
      ELSE
        tmp=t*MeshVel(SlidingMeshDirection(iPartition),iPartition)/delta_phi
      END IF
    CASE DEFAULT
      CALL Abort(__STAMP__, &
                       'Unknown SlidingMeshType in PrepareSM.')
  END SELECT

  SMo(iInterface)=tmp-REAL(INT(tmp))                             !fraction of current side in rotation
  SMn(iInterface)=MOD(INT(tmp),nAzimuthalSides(iGlobInterface))  !number of completed sides in rotation

  ! If we move in negative direction, correct SMo and SMn by an SM period
  IF(SMo(iInterface).LT.0) THEN
    SMo(iInterface) = SMo(iInterface)+1.
    SMn(iInterface) = SMn(iInterface)-1+nAzimuthalSides(iGlobInterface)
  END IF

  ! UPDATE PROJECTION / INTERPOLATION MATRICES
  !----------------------------------------------------------------------------------------------------
  CALL MortarBasis_BigToSmall(0,PP_N,NodeType, M_SM_0_1(:,:,iInterface), M_SM_0_2(:,:,iInterface), SMo_In=SMo(iInterface))
  CALL MortarBasis_SmallToBig(0,PP_N,NodeType, M_SM_1_0(:,:,iInterface), M_SM_2_0(:,:,iInterface), SMo_In=SMo(iInterface))

  M_SM_1_0(:,:,iInterface)=M_SM_1_0(:,:,iInterface)*    SMo(iInterface)  !weighting of projection matrices with lentgh
  M_SM_2_0(:,:,iInterface)=M_SM_2_0(:,:,iInterface)*(1.-SMo(iInterface)) !weighting of projection matrices with lentgh
END DO ! iInterface


! UPDATE NEIGHBOUR INFORMATION:
!----------------------------------------------------------------------------------------------------
IF(ANY(SMn.NE.SMn_prev))THEN !neighbour information has changed
  SM_nMPIMortars_Proc(:)=0

  DO i=STAT,ROT
  IF((i.EQ.STAT).AND.(.NOT.IAmAStatProc)) CYCLE
  IF((i.EQ.ROT ).AND.(.NOT.IAmARotProc )) CYCLE
    DO iSide=firstSMSide,lastSMSide
      DO iInterface=1,nSlidingMeshInterfaces
        iGlobInterface = LocToGlobInterface(iInterface)
        IF(SlidingMeshInfo(5,iSide).NE.iGlobInterface) CYCLE
        DO iMortar=1,2
          !get globally unique indices
          iiMortar=2*(iSide-firstSMSide)+iMortar
          iL = SlidingMeshInfo(2,iSide)
          iAz= SlidingMeshInfo(1,iSide)
          IF(i.EQ.ROT)THEN
            iAz=iAz+SMn(iInterface)+(2-iMortar)
            IF(iAz.GT.nAzimuthalSides(iGlobInterface)) iAz=iAz-nAzimuthalSides(iGlobInterface)
            iProc=AL_ToStatProc(iGlobInterface,iAz,iL)
          ELSE
            iAz=iAz-SMn(iInterface)-(2-iMortar)
            IF(iAz.LT.1) iAz=iAz+nAzimuthalSides(iGlobInterface)
            iProc=AL_ToRotProc(iGlobInterface,iAz,iL)
            iAz=iAz+SMn(iInterface)+(2-iMortar)
            IF(iAz.GT.nAzimuthalSides(iGlobInterface)) iAz=iAz-nAzimuthalSides(iGlobInterface)
          END IF
          !get a matrix with mapping to globally unique indices
          SideToSMMortar(iiMortar,S2SM_IPROC     )=iProc
          SideToSMMortar(iiMortar,S2SM_IINTERFACE)=iGlobInterface
          SideToSMMortar(iiMortar,S2SM_IAZ       )=iAz
          SideToSMMortar(iiMortar,S2SM_IL        )=iL
          SideToSMMortar(iiMortar,S2SM_IMORTAR   )=iMortar
          SideToSMMortar(iiMortar,S2SM_ISIDE     )=iSide
          SM_nMPIMortars_Proc(iProc)=SM_nMPIMortars_Proc(iProc)+1
        END DO !iMortar
      END DO !iInterface
    END DO !iSide

    !sort by proc and by globally unique inidces:
    !first by proc, then by interface, then by azimuth, then by layer, then by iMortar.
    !this allows the MPI receiving proc to know the sorting a priori
    CALL QSort4int(SideToSMMortar)

    !inverse mapping (which is eventually used)
    DO iiMortar=1,2*nSMSides
      iMortar = SideToSMMortar(iiMortar,S2SM_IMORTAR)
      iSide   = SideToSMMortar(iiMortar,S2SM_ISIDE  )
      IF(i.EQ.STAT)THEN
        SMMortarToStatSide(iMortar,iSide) = iiMortar
      ELSE
        SMMortarToRotSide(iMortar,iSide) = iiMortar
      END IF
    END DO !iiMortar
  END DO !i (STAT/ROT)
END IF !(SMn.NE.SMn_prev)
SMn_prev=SMn


! UPDATE NORMAL/TANGENT VECTORS
!----------------------------------------------------------------------------------------------------
IF (IAmAStatProc) THEN
  DO iInterface=1,nSlidingMeshInterfaces
    M1=>M_SM_0_1(:,:,iInterface); M2=>M_SM_0_2(:,:,iInterface)
    DO iSide=firstSMSide,lastSMSide

      ! Cycle if current SM side is not on iInterface
      IF (SlidingMeshInfo(SM_SIDE_TO_INTERFACE,iSide).NE.LocToGlobInterface(iInterface)) CYCLE

      Mortars(:) = SMMortarToStatSide(:,iSide)

      ! Interpolate stator side's normal/tangent vectors onto sm mortars
      CALL DoFlip3(NormVec(:,:,:,0,iSide),VecTmp,SMFlip(iSide))
      CALL InterpolateSM3(VecTmp,NormVec_SM( :,:,:,Mortars(1)),NormVec_SM( :,:,:,Mortars(2)),M1,M2,isRotor=.FALSE.)
      CALL DoFlip3(TangVec1(:,:,:,0,iSide),VecTmp,SMFlip(iSide))
      CALL InterpolateSM3(VecTmp,TangVec1_SM(:,:,:,Mortars(1)),TangVec1_SM(:,:,:,Mortars(2)),M1,M2,isRotor=.FALSE.)
      CALL DoFlip3(TangVec2(:,:,:,0,iSide),VecTmp,SMFlip(iSide))
      CALL InterpolateSM3(VecTmp,TangVec2_SM(:,:,:,Mortars(1)),TangVec2_SM(:,:,:,Mortars(2)),M1,M2,isRotor=.FALSE.)

    END DO ! iSide
  END DO ! iInterface
END IF

END SUBROUTINE PrepareSM



!============================================================================================================================
!> Determines the global flip of a side (i.e. the p-q-orientation of the master system in global x-y-z-coordinates)
!>
!>             |            |            |            |            |            |            |            |            |
!>   CoordSys  |   flip=1   |   flip=2   |   flip=3   |   flip=4   |   flip=5   |   flip=6   |   flip=7   |   flip=8   |
!>             |            |            |            |            |            |            |            |            |
!>  Az ^       |  p ^       |       ^ p  |        q   |   q        |  q ^       |       ^ q  |        p   |   p        |
!>     |       |    |       |       |    |    +--->   |   <---+    |    |       |       |    |    +--->   |   <---+    |
!>     +--->   |    +--->   |   <---+    |    |       |       |    |    +--->   |   <---+    |    |       |       |    |
!>         L   |        q   |   q        |  p v       |       v p  |        p   |   p        |  q v       |       v q  |
!>             |            |            |            |            |            |            |            |            |
!>
!============================================================================================================================
SUBROUTINE GetSMFlip(SideID,flip_out)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Mesh_Vars,           ONLY: Face_xGP,IntToPart,SlidingMeshInfo
USE MOD_SM_Vars,             ONLY: SlidingMeshDirection,SlidingMeshType
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER,INTENT(IN ) :: SideID
INTEGER,INTENT(OUT) :: flip_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                :: xTmp(1:2,0:PP_N,0:PP_NZ)
REAL,PARAMETER      :: eps=1.E-9
INTEGER             :: i,j
INTEGER             :: iPartition,iSlidingMeshType
!============================================================================================================================
! Get sm partition of side
iPartition       = IntToPart(SlidingMeshInfo(SM_SIDE_TO_INTERFACE,SideID))
iSlidingMeshType = SlidingMeshType(iPartition)

SELECT CASE (iSlidingMeshType)
  CASE(SM_TYPE_ROTATION_RADIAL)
    xTmp(1,:,:)=ATAN2(-Face_xGP(2,:,:,0,SideID),-Face_xGP(1,:,:,0,SideID))      ! Azimuth
    xTmp(2,:,:)=Face_xGP(3,:,:,0,SideID)                                        ! z
  CASE(SM_TYPE_ROTATION_AXIAL)
    xTmp(1,:,:)=ATAN2(-Face_xGP(2,:,:,0,SideID),-Face_xGP(1,:,:,0,SideID))      ! Azimuth
    xTmp(2,:,:)=SQRT(Face_xGP(1,:,:,0,SideID)**2+Face_xGP(2,:,:,0,SideID)**2)   ! r
  CASE(SM_TYPE_PLANAR)
    xTmp(1,:,:)= Face_xGP(SlidingMeshDirection(iPartition),:,:,0,SideID) ! x or y
    xTmp(2,:,:)=-Face_xGP(3                               ,:,:,0,SideID) ! z
  CASE DEFAULT
    CALL Abort(__STAMP__, &
                    'Unknown SlidingMeshType in GetSMFlip.')
END SELECT


! Separate treatment for the side including ATAN2 discontinuity at (x>0.,y=0.)
! ATTENTION: only valid for elements covering not more than 180° of sm interface!
IF ((iSlidingMeshType.EQ.SM_TYPE_ROTATION_RADIAL).OR.(iSlidingMeshType.EQ.SM_TYPE_ROTATION_AXIAL)) THEN
  ! Check whether side crosses x-axis (sign change in y) for x>0.
  IF (     (ANY((Face_xGP(2,:,:,0,SideID).LT.0.).AND.(Face_xGP(1,:,:,0,SideID).GT.0.)))  &
      .AND.(ANY((Face_xGP(2,:,:,0,SideID).GE.0.).AND.(Face_xGP(1,:,:,0,SideID).GT.0.)))) THEN
    DO i=0,PP_NZ;DO j=0,PP_N
      ! Add 2*PI to counteract the jump from +PI to -PI within side
      IF (xTmp(1,j,i).LT.0.) xTmp(1,j,i)=xTmp(1,j,i)+2*PP_PI
    END DO; END DO ! i,j
  END IF
END IF


#if PP_dim==3
IF      ( (xTmp(1,PP_N,0)-xTmp(1,0,0).GT. eps) .AND. (ABS(xTmp(1,0,PP_N)-xTmp(1,0,0)).LT.eps) ) THEN
  IF    ( (xTmp(2,0,PP_N)-xTmp(2,0,0).GT. eps) .AND. (ABS(xTmp(2,PP_N,0)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=1
  ELSEIF( (xTmp(2,0,PP_N)-xTmp(2,0,0).LT.-eps) .AND. (ABS(xTmp(2,PP_N,0)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=2
  ELSE
    CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
  END IF
ELSEIF  ( (xTmp(1,PP_N,0)-xTmp(1,0,0).LT.-eps) .AND. (ABS(xTmp(1,0,PP_N)-xTmp(1,0,0)).LT.eps) ) THEN
  IF    ( (xTmp(2,0,PP_N)-xTmp(2,0,0).GT. eps) .AND. (ABS(xTmp(2,PP_N,0)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=3
  ELSEIF( (xTmp(2,0,PP_N)-xTmp(2,0,0).LT.-eps) .AND. (ABS(xTmp(2,PP_N,0)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=4
  ELSE
    CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
  END IF
ELSEIF  ( (xTmp(1,0,PP_N)-xTmp(1,0,0).GT. eps) .AND. (ABS(xTmp(1,PP_N,0)-xTmp(1,0,0)).LT.eps) ) THEN
  IF    ( (xTmp(2,PP_N,0)-xTmp(2,0,0).GT. eps) .AND. (ABS(xTmp(2,0,PP_N)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=5
  ELSEIF( (xTmp(2,PP_N,0)-xTmp(2,0,0).LT.-eps) .AND. (ABS(xTmp(2,0,PP_N)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=6
  ELSE
    CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
  END IF
ELSEIF  ( (xTmp(1,0,PP_N)-xTmp(1,0,0).LT.-eps) .AND. (ABS(xTmp(1,PP_N,0)-xTmp(1,0,0)).LT.eps) ) THEN
  IF    ( (xTmp(2,PP_N,0)-xTmp(2,0,0).GT. eps) .AND. (ABS(xTmp(2,0,PP_N)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=7
  ELSEIF( (xTmp(2,PP_N,0)-xTmp(2,0,0).LT.-eps) .AND. (ABS(xTmp(2,0,PP_N)-xTmp(2,0,0)).LT.eps) ) THEN
    flip_out=8
  ELSE
    CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
  END IF
ELSE
  CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
END IF
#else
IF     (xTmp(1,PP_N,0)-xTmp(1,0,0).GT. eps) THEN
  flip_out=1
ELSEIF (xTmp(1,PP_N,0)-xTmp(1,0,0).LT.-eps) THEN
  flip_out=3
ELSE
  CALL Abort(__STAMP__,'Error: could not determine SM Side Flip' )
END IF
#endif /*PP_dim==3*/
END SUBROUTINE GetSMFlip



!===================================================================================================================================
!> QSort4int:
!> Uses the Quicksort algorithm to sort an INTEGER table with five relevant columns.
!> There may be an arbitrary number of additional columns that will be moved around at the same
!> time but do not influence the sorting process (except for slowing things down of course!)
!===================================================================================================================================
PPURE RECURSIVE SUBROUTINE Qsort4int(A)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: iq  ! ?
!===================================================================================================================================
IF(size(A,1) .GT. 1) THEN
   CALL Partition4int(A, iq)
   CALL Qsort4int(A(:iq-1,:))
   CALL Qsort4int(A(iq:,:))
END IF
END SUBROUTINE Qsort4int



!===================================================================================================================================
!>  Sorting routine used by QSort4int above. This routine is PRIVATE
!===================================================================================================================================
PPURE SUBROUTINE Partition4int(A, marker)
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT/ OUTPUT VARIABLES
INTEGER, INTENT(INOUT), DIMENSION(:,:) :: A   ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
INTEGER, INTENT(OUT)                   :: marker  ! ?
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: i,j                ! ?
INTEGER   :: temp(size(A,2))        ! ?
INTEGER   :: x(5)              ! pivot point
!===================================================================================================================================
x(:) = A(1,1:5)
i= 0
j= size(A,1) + 1

DO
  j = j-1
  DO
    IF(A(j,1) < x(1))THEN
      EXIT
    ELSE IF(A(j,1) == x(1))THEN
      IF(A(j,2) < x(2))THEN
        EXIT
      ELSE IF(A(j,2) == x(2))THEN
        IF(A(j,3) < x(3))THEN
          EXIT
        ELSE IF(A(j,3) == x(3))THEN
          IF(A(j,4) < x(4))THEN
            EXIT
          ELSEIF(A(j,4) == x(4))THEN
            IF(A(j,5) <= x(5))THEN
              EXIT
            END IF
          END IF
        END IF
      END IF
    END IF
    j = j-1
  END DO
  i = i+1
  DO
    IF(A(i,1) > x(1))THEN
      EXIT
    ELSE IF(A(i,1) == x(1))THEN
      IF(A(i,2) > x(2))THEN
        EXIT
      ELSE IF(A(i,2) == x(2))THEN
        IF(A(i,3) > x(3))THEN
          EXIT
        ELSE IF(A(i,3) == x(3))THEN
          IF(A(i,4) > x(4))THEN
            EXIT
          ELSE IF(A(i,4) == x(4))THEN
            IF(A(i,5) >= x(5))THEN
              EXIT
            END IF
          END IF
        END IF
      END IF
    END IF
    i = i+1
  END DO
  IF (i < j) THEN
     ! exchange A(i) and A(j)
     temp(:) = A(i,:)
     A(i,:) = A(j,:)
     A(j,:) = temp(:)
  ELSE IF (i == j) THEN
     marker = i+1
     RETURN
  ELSE
     marker = i
     RETURN
  END IF
END DO
END SUBROUTINE Partition4int



!==================================================================================================================================
!> Deallocate SM variables
!==================================================================================================================================
SUBROUTINE FinalizeSM()
! MODULES
USE MOD_SM_Vars
IMPLICIT NONE
!==================================================================================================================================
SMInitIsDone=.FALSE.

SDEALLOCATE(SlidingMeshType)
SDEALLOCATE(SlidingMeshDirection)
SDEALLOCATE(SlidingMeshBoundaries)
SDEALLOCATE(SMn_prev)
SDEALLOCATE(M_SM_0_1)
SDEALLOCATE(M_SM_0_2)
SDEALLOCATE(M_SM_1_0)
SDEALLOCATE(M_SM_2_0)
SDEALLOCATE(NormVec_SM )
SDEALLOCATE(TangVec1_SM)
SDEALLOCATE(TangVec2_SM)
SDEALLOCATE(U_MorStat)
SDEALLOCATE(U_MorRot)
SDEALLOCATE(UPrim_MorStat)
SDEALLOCATE(UPrim_MorRot )
SDEALLOCATE(Flux_MorStat)
SDEALLOCATE(Flux_MorRot)
SDEALLOCATE(SMFlip)
SDEALLOCATE(SMMortarToStatSide)
SDEALLOCATE(SMMortarToRotSide)
SDEALLOCATE(SM_nMPIMortars_Proc)
SDEALLOCATE(Face_vGPSM)
#if PARABOLIC
SDEALLOCATE(gradUx_MorStat)
SDEALLOCATE(gradUy_MorStat)
SDEALLOCATE(gradUz_MorStat)
SDEALLOCATE(gradUx_MorRot )
SDEALLOCATE(gradUy_MorRot )
SDEALLOCATE(gradUz_MorRot )
#if EDDYVISCOSITY
SDEALLOCATE(muSGS_MorStat)
SDEALLOCATE(muSGS_MorRot )
#endif /*EDDYVISCOSITY*/
#endif /*PARABOLIC*/
END SUBROUTINE FinalizeSM

END MODULE MOD_SM










