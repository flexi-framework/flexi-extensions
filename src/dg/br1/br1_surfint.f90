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
#if PARABOLIC
#include "flexi.h"

!==================================================================================================================================
!> \brief Contains the surface integral routine for the BR1 lifting operation.
!>
!> This module contains the subroutine that performs the surface integral operation for the BR1 lifting. The BR1 lifting is
!> impelemented in strong or weak form. The fluxes for the strong form are different for the slave and the master sides.
!==================================================================================================================================
MODULE MOD_Lifting_SurfInt
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE Lifting_SurfInt
  MODULE PROCEDURE Lifting_SurfInt
END INTERFACE


PUBLIC::Lifting_SurfInt
!==================================================================================================================================
CONTAINS

!==================================================================================================================================
!> \brief Surface integral in the BR1 scheme optimized for performance, for weak or strong formulation.
!>
!> Performs the surface integral for the BR1 routine. Uses the DoSurfInt routine from the DG operator to perform actual
!> integration.
!> If we use the strong formulation, the inner solution is substracted from the numerical flux. This is done in the FillFlux
!> routines for the master side. The flux on the master side is \f$ \frac{1}{2} (U^+ + U^-) - U^- = \frac{1}{2} (U^+ - U^-) \f$
!> since \f$ \frac{1}{2} (U^+ + U^-) \f$ is the numerical flux in the BR1 scheme. On the slave side, the flux becomes
!> \f$ \frac{1}{2} (U^+ + U^-) - U^+ = \frac{1}{2} (-U^+ + U^-) \f$ which is simply the master side flux multiplied by
!> \f$ -1 \f$. This means we don't have to flip the sign on the flux for the slave side in strong form as we normally do to
!> get the flux on the slave side.
!==================================================================================================================================
SUBROUTINE Lifting_SurfInt(Nloc,Flux_master,Flux_slave,gradU,doMPISides,L_HatMinus,L_HatPlus,weak)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_SurfintPrim,        ONLY: DoSurfIntPrim
USE MOD_Mesh_Vars,          ONLY: SideToElem,nSides,nElems
USE MOD_Mesh_Vars,          ONLY: firstMPISide_YOUR,lastMPISide_MINE
USE MOD_Mesh_Vars,          ONLY: firstSMSide,lastSMSide
USE MOD_Mesh_Vars,          ONLY: S2V2
USE MOD_Mesh_Vars,          ONLY: nElems
#if FV_ENABLED
USE MOD_FV_Vars,            ONLY: FV_Elems_master,FV_Elems_slave
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)     :: Nloc                                             !< polynomial degree
LOGICAL,INTENT(IN)     :: doMPISides                                       !< = .TRUE. only MPISides_YOUR+MPIMortar are filled
                                                                           !< =.FALSE. BCSides+(Mortar-)InnerSides+MPISides_MINE
REAL,TARGET,INTENT(IN) :: Flux_master(1:PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides)   !< lifting flux (identical for master and slace)
REAL,TARGET,INTENT(IN) :: Flux_slave( 1:PP_nVarPrim,0:Nloc,0:ZDIM(Nloc),nSides)   !< slave lifting flux, only for sm without MPI
                                                                                  !< Otherwise dummy argument
REAL,INTENT(IN)        :: L_HatPlus(0:Nloc)                                !< lagrange polynomials at xi=+1 and pre-divided by
                                                                           !< integration weight
REAL,INTENT(IN)        :: L_HatMinus(0:Nloc)                               !< lagrange polynomials at xi=-1 and pre-divided by
                                                                           !< integration weight
REAL,INTENT(INOUT)     :: gradU(PP_nVarPrim,0:Nloc,0:Nloc,0:ZDIM(Nloc),1:nElems) !< time derivative of solution
LOGICAL,INTENT(IN)     :: weak                                             !< switch for weak or strong formulation
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: ElemID,nbElemID,locSideID,nblocSideID,SideID,p,q,flip
INTEGER            :: firstSideID,lastSideID,firstCycleID,lastCycleID
REAL               :: FluxTmp(1:PP_nVarPrim,0:Nloc,0:ZDIM(Nloc))
REAL,POINTER       :: Flux(:,:,:,:)
!==================================================================================================================================
IF(doMPISides)THEN       ! Only start with SM sides if there are any.
     firstSideID = MERGE(firstSMSide,firstMPISide_YOUR,lastSMSide.GE.firstSMSide)
      lastSideID = nSides
    firstCycleID = lastSMSide+1
     lastCycleID = firstMPISide_YOUR-1
ELSE
     firstSideID = 1
      lastSideID = lastMPISide_MINE
#if USE_MPI
    firstCycleID = firstSMSide
     lastCycleID = lastSMSide
#endif
END IF

! Generally use master flux, as master and slave lifting fluxes are identical (not on sm sides!)
Flux => Flux_master

DO SideID=firstSideID,lastSideID
#if USE_MPI
  IF((SideID.GE.firstCycleID).AND.(SideID.LE.lastCycleID)) CYCLE
#endif

  ElemID      = SideToElem(S2E_ELEM_ID,   SideID)
  nbElemID    = SideToElem(S2E_NB_ELEM_ID,SideID)

  ! master sides
  IF(ElemID.GT.0)THEN
    IF (FV_Elems_master(SideID).EQ.0) THEN ! DG element

#if !(USE_MPI)
      ! Reset flux pointer to master flux
      Flux => Flux_master
#endif

      locSideID = SideToElem(S2E_LOC_SIDE_ID,SideID)
      flip      = 0
      ! orient flux to fit flip and locSide to element local system
      DO q=0,ZDIM(Nloc); DO p=0,Nloc
        ! note: for master sides, the mapping S2V2 should be a unit matrix
        FluxTmp(:,S2V2(1,p,q,flip,locSideID),S2V2(2,p,q,flip,locSideID)) = Flux(:,p,q,SideID)
      END DO; END DO ! p,q
#if   (PP_NodeType==1)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      locSideID,gradU(:,:,:,:,ElemID))
#elif (PP_NodeType==2)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),locSideID,gradU(:,:,:,:,ElemID))
#endif
    END IF
  END IF

  ! slave sides
  IF(nbElemID.GT.0)THEN
    IF (FV_Elems_slave(SideID).EQ.0) THEN ! DG element

#if !(USE_MPI)
      ! Without MPI, the proc has the master and slave SM sides, which have different lifting
      ! fluxes, due to their changing position over time.  We therefore have to differentiate
      ! between master and slave sm lifting fluxes for use without MPI.
      IF ((SideID.GE.firstSMSide).AND.(SideID.LE.lastSMSide)) THEN
        Flux => Flux_slave
      END IF
#endif

      nblocSideID = SideToElem(S2E_NB_LOC_SIDE_ID,SideID)
      flip        = SideToElem(S2E_FLIP,SideID)
      ! orient flux to fit flip and locSide to element local system
      IF(weak)THEN
        DO q=0,ZDIM(Nloc); DO p=0,Nloc
          ! p,q are in the master RHS system, they need to be transformed to the slave volume system using S2V2 mapping
          FluxTmp(:,S2V2(1,p,q,flip,nblocSideID),S2V2(2,p,q,flip,nblocSideID))=-Flux(:,p,q,SideID)
        END DO; END DO ! p,q
      ELSE
        ! In strong form, don't flip the sign since the slave flux is the negative of the master flux
        DO q=0,ZDIM(Nloc); DO p=0,Nloc
          ! p,q are in the master RHS system, they need to be transformed to the slave volume system using S2V2 mapping
          FluxTmp(:,S2V2(1,p,q,flip,nblocSideID),S2V2(2,p,q,flip,nblocSideID))= Flux(:,p,q,SideID)
        END DO; END DO ! p,q
      END IF
#if   (PP_NodeType==1)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus,   L_HatPlus,      nblocSideID,gradU(:,:,:,:,nbElemID))
#elif (PP_NodeType==2)
      CALL DoSurfIntPrim(Nloc,FluxTmp,L_HatMinus(0),L_HatPlus(Nloc),nblocSideID,gradU(:,:,:,:,nbElemID))
#endif
    END IF
  END IF
END DO ! SideID=1,nSides
END SUBROUTINE Lifting_SurfInt

END MODULE MOD_Lifting_SurfInt
#endif /*PARABOLIC*/
