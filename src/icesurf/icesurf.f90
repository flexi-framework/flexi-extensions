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
#include "eos.h"

!==================================================================================================================================
!> Send data to ice surface
!==================================================================================================================================
MODULE MOD_IceSurf
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitIceSurf
  MODULE PROCEDURE InitIceSurf
END INTERFACE

INTERFACE CalcIceSurfData
  MODULE PROCEDURE CalcIceSurfData
END INTERFACE

INTERFACE FinalizeIceSurf
  MODULE PROCEDURE FinalizeIceSurf
END INTERFACE

PUBLIC :: DefineParametersIceSurf
PUBLIC :: InitIceSurf
PUBLIC :: CalcIceSurfData
PUBLIC :: FinalizeIceSurf
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters
!==================================================================================================================================
SUBROUTINE DefineParametersIceSurf()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("IceSurf")
CALL prms%CreateLogicalOption('doCalcIceSurfData'  , "Calculate and write ice surface data for icing model","T")
CALL prms%CreateIntOption('NOutSurf'  , "Surface FV Output polynomial degree")
END SUBROUTINE DefineParametersIceSurf


!=================================================================================================================================
!> Initialize ice surf.
!=================================================================================================================================
SUBROUTINE InitIceSurf()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_IceSurf_Vars
USE MOD_ReadInTools,        ONLY:GETLOGICAL,GETINT
USE MOD_Mesh_Vars, ONLY: nBCSides,BoundaryType,BC,ElemToSide,nElems
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
!LOCAL VARIABLES
INTEGER             :: SideID, BCType, iElem, LocSideID
INTEGER             :: nWallSidesProcs(0:nProcessors-1)
!==================================================================================================================================
IF (IceSurfInitIsDone) THEN
  CALL CollectiveStop(__STAMP__,&
    'InitIceSurf already called.')
END IF
SWRITE(UNIT_StdOut,'(132("-"))')
SWRITE(UNIT_stdOut,'(A)') ' INIT ICESURF...'

doCalcIceSurfData = GETLOGICAL("doCalcIceSurfData","T")

IF(doCalcIceSurfData)THEN 
  nWallSides = 0
  DO SideID=1,nBCSides
    BCType  = BoundaryType(BC(SideID),BC_TYPE)
    IF(ANY(BCType.EQ.WALLBCTYPES())) nWallSides = nWallSides +1 
  END DO

#if USE_MPI
  CALL MPI_ALLGATHER(nWallSides,1,MPI_INTEGER,nWallSidesProcs,1,MPI_INTEGER,MPI_COMM_FLEXI,iError)
  nWallSidesGlob = SUM(nWallSidesProcs)
  offsetWallSides = SUM(nWallSidesProcs(0:myRank-1))
#else /*USE_MPI*/
  nWallSidesGlob = nWallSides
  offsetWallSides = 0
#endif /*USE_MPI*/

  !Map to sort by ElemID, which is global, whereas SideID is process-local.
  ALLOCATE(MapIceSurf(nWallSides))
  nWallSides = 0
  DO iElem=1,nElems
#if PP_dim == 3
    DO LocSideID=1,6
#else
    DO LocSideID=2,5
#endif
      SideID = ElemToSide(E2S_SIDE_ID,LocSideID,iElem)
      IF(SideID.GT.nBCSides) CYCLE
      BCType  = Boundarytype(BC(SideID),BC_TYPE)
      IF(.NOT.ANY(BCType.EQ.WALLBCTYPES())) CYCLE
      nWallSides = nWallSides + 1 
      MapIceSurf(nWallSides) = SideID
    END DO 
  END DO 

  NOutSurf = GETINT("NOutSurf")
  ALLOCATE(IceSurfData(ICS_NVAR,0:NOutSurf,0:ZDIM(NOutSurf),nWallSides))
ENDIF 


IceSurfInitIsDone = .TRUE.
SWRITE(UNIT_stdOut,'(A)')' INIT ICESURF DONE!'
SWRITE(UNIT_StdOut,'(132("-"))')
END SUBROUTINE InitIceSurf


!==================================================================================================================================
!> Note: For truly instationary results, this should be called after the first DGTimeDerivative of the time step to have 
!>       consistent U and Ut
!==================================================================================================================================
SUBROUTINE CalcIceSurfData()
! MODULES
USE MOD_PreProc
USE MOD_IceSurf_Vars
USE MOD_DG_Vars,             ONLY:UPrim_master
USE MOD_Mesh_Vars,           ONLY:Face_xGP
USE MOD_Interpolation_Vars,  ONLY:NodeType
USE MOD_ChangeBasisByDim,    ONLY:ChangeBasisSurf
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: SideID, WallSideID
REAL                           :: Vdm(0:NOutSurf,0:PP_N)   !< Vandermonde matrix for converstion from DG to FV
REAL                           :: IceSurfDG(ICS_NVAR,0:PP_N,0:PP_NZ)
!==================================================================================================================================
CALL ICE_FV_GetVandermonde(PP_N,NOutSurf,NodeType,Vdm)


DO WallSideID=1,nWallSides
  SideID = MapIceSurf(WallSideID)
  !TODO: This is only correct for 1D surf!!
  !IceSurfDG(ICS_VCEL,:,:) = SurfElem(:,:,0,SideID) *2. / (PP_N+1)

  IceSurfDG(ICS_X   ,:,:) = Face_xGP(1,:,:,0,SideID)
  IceSurfDG(ICS_Y   ,:,:) = Face_xGP(2,:,:,0,SideID)
  IceSurfDG(ICS_P,:,:) = UPrim_master(5,:,:,SideID)

  !TODO: Only change basis if not FVElem
  CALL ChangeBasisSurf(ICS_NVAR,PP_N,NOutSurf,Vdm,IceSurfDG,IceSurfData(:,:,:,WallSideID))
END DO

!==================================================================================================================================
! TODO: 
! fix bug where first entry of IceSurfData(ICS_X...) is 0
! rename all "ice" in varnames 
!==================================================================================================================================

END SUBROUTINE CalcIceSurfData



!==================================================================================================================================
!> COPIED FROM FV_BASIS
!> Build Vandermondes to switch solution between DG and FV sub-cells.
!==================================================================================================================================
SUBROUTINE ICE_FV_GetVandermonde(N_in,N_Out,NodeType_in,FV_Vdm)
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_Interpolation ,ONLY: GetVandermonde,GetNodesAndWeights
USE MOD_Basis         ,ONLY: InitializeVandermonde
USE MOD_Mathtools     ,ONLY: INVERSE
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)            :: N_in                    !< Number of 1D input points / output points
INTEGER,INTENT(IN)            :: N_Out                   !< Number of 1D input points / output points
CHARACTER(LEN=255),INTENT(IN) :: NodeType_in             !< Type of 1D input points
REAL,INTENT(OUT)              :: FV_Vdm(0:N_Out,0:N_in)   !< Vandermonde matrix for converstion from DG to FV
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                   :: FV_X(0:N_Out),FV_w,FV_BdryX(0:N_Out+1)
REAL,DIMENSION(0:N_In) :: xGP,wGP,wBary
REAL                   :: SubxGP(0:N_In)
REAL                   :: VDM(0:N_In,0:N_In)
INTEGER                :: i,j,k
!==================================================================================================================================
CALL GetNodesAndWeights(N_in,NodeType_in,xGP,wGP,wBary)

! one DG cell [-1,1] is divided in FV-Subcells
!
! -1                  i-th Subcell                              1
!  |------- ... ----|---x---x---x---|---- ... -----|------------|
!                     x = xGP in Subcell
!
!  xFV_i(k) = k-th Gauss point in i-th Subcell
!
! We have to integrate the solution U in each Subcell (loop over i) and compute with this the integral mean value, which
! then ist the Finite-Volume solution U_FV of this Subcell.
! Therefore we must compute the solution in each Gauss point xFV_i of the Subcell, which is done by evaluating all Lagrange
! polynomials l_j in all Gauss points xFV_i => stored in the matrix VDM_i = [ l_j(xFV_i(k)) ]_kj
! Multiplying this Vandermonde with the DG-solution U gives us the solution in the Gauss points xFV_i of the i-th Subcell:
!    u_i = VDM_i . U         meaning u_i(k) = U(xFV_i(k))
! The integration over the Subcell is done by Gauss-Quadrature:
!    \int_{Subcell(i)} U dx = [ \sum_{k=0..N} u_i(k) . wGP(k) ] * a       with a = (width of Subcell(i)) / 2
! ( The /2 comes from the width of the reference element. )
! We get the integral mean value and therewith the Finite-Volume solution of the i-th Subcell  by dividing this integral by
! the width of the i-th Subcell :
!    u_FV(i) = [ \sum_{k=0..N} wGP(k) * u_i(k) ] / 2
!
! All this can be write in Matrix-Vector-Multiplication:
!    u_FV(i) = 1/2 * wGP^T . VDM_i . U
! If we combine the vectors (wGP^T . VDM_i) for the i-th Subcell for all Subcells in one matrix we get :
!    u_FV = 1/2 * FV_Vdm . U          with FV_Vdm = ( wGP^T.VDM_0 // wGP^T.VDM_1 // ... // wGP^T. VDM_N )
! With the inverse of the matrix FV_Vdm we then can reconstruct the DG-Solution from a FV-Solution
!
! Algorithm :
FV_Vdm = 0.0
CALL ICE_FV_Build_X_w_BdryX(N_Out, FV_X, FV_w, FV_BdryX)
!FV_BdryX = (FV_BdryX + 1.)/2. -1.
DO i=0,N_Out
  ! 1. Compute the Gauss points xFV_i in the i-th Subcell
  SubxGP = FV_BdryX(i) + (xGP + 1.)/2. * (FV_BdryX(i+1) - FV_BdryX(i))
  ! 2. Evaluate the all Lagrange-Polys in all Gauss points xFV_i of the i-th Subcell  =>  store in VDM
  CALL InitializeVandermonde(N_in,N_in,wBary,xGP,SubxGP,VDM(:,:))
  ! 3. Multiply wGP^T with VDM and store it in the i-th row of the matrix FV_Vdm
  DO j=0,N_in
    DO k=0,N_in
      FV_Vdm(i,j) = FV_Vdm(i,j) + wGP(k) * VDM(k,j)
    END DO
  END DO
END DO
! 4. don't forget the 1/2
FV_Vdm = FV_Vdm * 0.5

END SUBROUTINE ICE_FV_GetVandermonde

!==================================================================================================================================
!> COPIED FROM FV_BASIS
!> Build positions FV_X, widths, and boundary positions
!==================================================================================================================================
! MODULES
SUBROUTINE ICE_FV_Build_X_w_BdryX(N_FV, FV_X, FV_w, FV_BdryX)
USE MOD_Basis
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_FV               !< polynomial degree of DG elements / number of sub-cells per direction (N+1)
REAL,INTENT(OUT)   :: FV_X(0:N_FV)       !< cell-centers of the sub-cells in reference space
REAL,INTENT(OUT)   :: FV_w            !< width of the sub-cells in reference space
REAL,INTENT(OUT)   :: FV_BdryX(0:N_FV+1) !< positions of the boundaries of the sub-cells in reference space
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: i
!==================================================================================================================================
FV_w  = 2.0 / (N_FV+1) ! equidistant widths of FV-Subcells

! calculate boundaries and nodes (midpoints) of FV-Subcells
FV_BdryX(0) = -1.0
DO i=1,N_FV+1
  FV_BdryX(i) = FV_BdryX(i-1) + FV_w
  FV_X(i-1)   = (FV_BdryX(i-1) + FV_BdryX(i))/2.0
END DO
END SUBROUTINE ICE_FV_Build_X_w_BdryX


!============================================================================================================================
!> Deallocate all global ice surf variables
!============================================================================================================================
SUBROUTINE FinalizeIceSurf()
! MODULES
USE MOD_IceSurf_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
!input parameters
!----------------------------------------------------------------------------------------------------------------------------
!output parameters
!----------------------------------------------------------------------------------------------------------------------------
!local variables
!============================================================================================================================
! Deallocate global variables
SDEALLOCATE(IceSurfData)
SDEALLOCATE(MapIceSurf)

IceSurfInitIsDone = .FALSE.
END SUBROUTINE FinalizeIceSurf


END MODULE MOD_IceSurf
