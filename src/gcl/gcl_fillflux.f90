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
!> Routine to fill the fluxes of the GCL.
!===================================================================================================================================
MODULE MOD_GCL_FillFlux
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------
INTERFACE GCL_FillFlux
  MODULE PROCEDURE GCL_FillFlux
END INTERFACE

PUBLIC::GCL_FillFlux
!===================================================================================================================================

CONTAINS


!===================================================================================================================================
!> Fill the fluxes of the GCL. The fluxes consist of the transformed negative mesh velocites, which are continuous 
!> across element faces. Because of this, no riemann solver is needed here.
!===================================================================================================================================
SUBROUTINE GCL_FillFlux(Flux_master,Flux_slave,doMPISides)
! MODULES
USE MOD_PreProc
USE MOD_Mesh_Vars,       ONLY: NormVec,SurfElem
USE MOD_Mesh_Vars,       ONLY: lastMPISide_MINE,firstMPISide_MINE,lastInnerSide,firstBCSide,nSides
USE MOD_Mesh_Vars,       ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_MoveMesh_Vars,   ONLY: Face_vGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
REAL,INTENT(OUT)   :: Flux_master(1,0:PP_N,0:PP_NZ,1:nSides)  !< GCL flux 
REAL,INTENT(OUT)   :: Flux_slave (1,0:PP_N,0:PP_NZ,1:nSides)  !< GCL flux
LOGICAL,INTENT(IN) :: doMPISides                             !< = .TRUE. only MINE MPISides are filled, =.FALSE. InnerSides  
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: SideID,p,q,firstSideID,lastSideID
!===================================================================================================================================
! fill flux for sides ranging between firstSideID and lastSideID, simply negative mesh velocity normal to interface 
IF(doMPISides)THEN
  ! fill only flux for MINE MPISides
  firstSideID = firstMPISide_MINE
  lastSideID  = lastMPISide_MINE
ELSE
  ! fill InnerSides and BCSides
  ! Inner Sides and BCSides can be treated in the same way in the GCL, since we need no boundary fluxes 
  firstSideID = firstBCSide
  lastSideID  = lastInnerSide
END IF

DO SideID=firstSideID,lastSideID
  ! Don't compute fluxes on big mortar sides (flux will be projected from small sides instead)
  IF ((SideID.GE.firstMortarInnerSide).AND.(SideID.LE.lastMortarInnerSide)) CYCLE
  DO q=0,PP_NZ
    DO p=0,PP_N
      Flux_master(1,p,q,SideID)=-(Face_vGP(1,p,q,SideID)*NormVec(1,p,q,0,SideID) + &
                                  Face_vGP(2,p,q,SideID)*NormVec(2,p,q,0,SideID) &
#if (PP_dim == 3)
                                 +Face_vGP(3,p,q,SideID)*NormVec(3,p,q,0,SideID) &
#endif
                                  )*SurfElem(p,q,0,SideID)
    END DO
  END DO
END DO ! SideID

! copy flux from master side to slave side
DO SideID=firstSideID ,lastSideID
  IF ((SideID.GE.firstMortarInnerSide).AND.(SideID.LE.lastMortarInnerSide)) CYCLE
  Flux_slave(:,:,:,SideID) = Flux_master(:,:,:,SideID)
END DO !SideID

END SUBROUTINE GCL_FillFlux

END MODULE MOD_GCL_FillFlux
