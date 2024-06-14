!=================================================================================================================================
! Copyright (c) 2010-2022  Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024  Prof. Andrea Beck
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://numericsresearchgroup.org
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
!> Contains routines to compute the riemann (Advection, Diffusion) for a given Face
!==================================================================================================================================
MODULE MOD_Riemann
! MODULES
IMPLICIT NONE
PRIVATE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  PURE SUBROUTINE RiemannInt(U_LL,U_RR,UB_LL,UB_RR,F)
    REAL,DIMENSION(PP_nVar),INTENT(IN) :: U_LL,U_RR
    REAL,DIMENSION(PP_nVarBase),INTENT(IN) :: UB_LL,UB_RR
    REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
  END SUBROUTINE
END INTERFACE

PROCEDURE(RiemannInt),POINTER :: Riemann_pointer    !< pointer defining the standard inner Riemann solver
PROCEDURE(RiemannInt),POINTER :: RiemannBC_pointer  !< pointer defining the standard BC    Riemann solver

INTEGER,PARAMETER      :: PRM_RIEMANN_SAME          = -1
INTEGER,PARAMETER      :: PRM_RIEMANN_LF            = 1
INTEGER,PARAMETER      :: PRM_RIEMANN_FLUXVECSPLIT  = 2
INTEGER,PARAMETER      :: PRM_RIEMANN_CENTRAL       = 3

INTERFACE InitRiemann
  MODULE PROCEDURE InitRiemann
END INTERFACE

INTERFACE Riemann
  MODULE PROCEDURE Riemann
END INTERFACE

INTERFACE FinalizeRiemann
  MODULE PROCEDURE FinalizeRiemann
END INTERFACE

PUBLIC::InitRiemann
PUBLIC::Riemann
PUBLIC::FinalizeRiemann
!==================================================================================================================================

PUBLIC::DefineParametersRiemann
CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES 
!==================================================================================================================================
CALL prms%SetSection("Riemann")
CALL prms%CreateIntFromStringOption('Riemann',   "Riemann solver to be used: LF, fluxvecsplit, central",&
                                                 "LF")
CALL addStrListEntry('Riemann','lf',             PRM_RIEMANN_LF)
CALL addStrListEntry('Riemann','fluxvecsplit',   PRM_RIEMANN_FLUXVECSPLIT)
CALL addStrListEntry('Riemann','central',        PRM_RIEMANN_CENTRAL)
CALL prms%CreateIntFromStringOption('RiemannBC', "Riemann solver used for boundary conditions: same, LF, fluxvecsplit",&
                                                 "same")
CALL addStrListEntry('RiemannBC','same',         PRM_RIEMANN_SAME)
CALL addStrListEntry('RiemannBC','lf',           PRM_RIEMANN_LF)
CALL addStrListEntry('RiemannBC','fluxvecsplit', PRM_RIEMANN_FLUXVECSPLIT)
CALL addStrListEntry('RiemannBC','central',      PRM_RIEMANN_CENTRAL)
END SUBROUTINE DefineParametersRiemann

!==================================================================================================================================!
!> Initialize Riemann solver routines, read inner and BC Riemann solver parameters and set pointers
!==================================================================================================================================!
SUBROUTINE InitRiemann()
! MODULES
USE MOD_Globals
USE MOD_ReadInTools ,ONLY: GETINTFROMSTR
!----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: Riemann
!==================================================================================================================================
Riemann = GETINTFROMSTR('Riemann')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_LF)
  Riemann_pointer => Riemann_LF
CASE(PRM_RIEMANN_FLUXVECSPLIT)
  Riemann_pointer => Riemann_FVS
CASE(PRM_RIEMANN_CENTRAL)
  Riemann_pointer => Riemann_central
!   CALL CollectiveStop(__STAMP__,&
!    'Flux vector splitting Riemann solver not yet implemented!')
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'Riemann solver not defined!')
END SELECT

Riemann = GETINTFROMSTR('RiemannBC')
SELECT CASE(Riemann)
CASE(PRM_RIEMANN_SAME)
  RiemannBC_pointer => Riemann_pointer
CASE(PRM_RIEMANN_LF)
  RiemannBC_pointer => Riemann_LF
CASE(PRM_RIEMANN_FLUXVECSPLIT)
  RiemannBC_pointer => Riemann_FVS
CASE(PRM_RIEMANN_CENTRAL)
  RiemannBC_pointer => Riemann_central
CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    'RiemannBC solver not defined!')
END SELECT

END SUBROUTINE InitRiemann

!==================================================================================================================================
!> Computes the numerical flux
!> Conservative States are rotated into normal direction in this routine and are NOT backrotatet: don't use it after this routine!!
!==================================================================================================================================
SUBROUTINE Riemann(Nloc,Fout,U_L,U_R,UB_L,UB_R,nv,t1,t2,doBC)
! MODULES
USE MOD_PreProc
USE MOD_Flux,   ONLY: EvalEulerFlux1D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                                    :: NLoc                    !< Polynomial degree
REAL,DIMENSION(PP_nVar,    0:Nloc,0:Nloc),INTENT(IN)  :: U_L                     !< Left state
REAL,DIMENSION(PP_nVar,    0:Nloc,0:Nloc),INTENT(IN)  :: U_R                     !< Right state
REAL,DIMENSION(PP_nVarBase,0:Nloc,0:Nloc),INTENT(IN)  :: UB_L                    !< baseflow state
REAL,DIMENSION(PP_nVarBase,0:Nloc,0:Nloc),INTENT(IN)  :: UB_R                    !< baseflow state
REAL,INTENT(IN)                                       :: nv(3,0:Nloc,0:Nloc)     !< Normal vector
REAL,INTENT(IN)                                       :: t1(3,0:Nloc,0:Nloc)     !< First tangential vector
REAL,INTENT(IN)                                       :: t2(3,0:Nloc,0:Nloc)     !< Second tangential vector
LOGICAL,INTENT(IN)                                    :: doBC                    !< Switch to do BC sides or not
REAL,INTENT(OUT)                                      :: Fout(PP_nVar,0:Nloc,0:Nloc)!< Flux
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                 :: i,j
REAL,DIMENSION(PP_nVar) :: U_LL,U_RR,F
REAL,DIMENSION(PP_nVarBase):: UB_LL,UB_RR
PROCEDURE(RiemannInt),POINTER :: Riemann_loc !< pointer defining the standard inner Riemann solver
!==================================================================================================================================
IF (doBC) THEN
  Riemann_loc => RiemannBC_pointer
ELSE
  Riemann_loc => Riemann_pointer
END IF

! Momentum has to be rotatet using the normal system individual for each
DO j=0,ZDIM(NLoc); DO i=0,NLoc
  ! left state: U_L
  U_LL(DENS)=U_L(DENS,i,j)
  ! rotate velocity in normal and tangential direction 
  U_LL(VEL1)=DOT_PRODUCT(U_L(VELV,i,j),nv(:,i,j))
  U_LL(VEL2)=DOT_PRODUCT(U_L(VELV,i,j),t1(:,i,j))
#if PP_dim==3
  U_LL(VEL3)=DOT_PRODUCT(U_L(VELV,i,j),t2(:,i,j))
#else
  U_LL(VEL3)=0.
#endif
  U_LL(PRES)=U_L(PRES,i,j)
  ! left baseflow state: UB_L
  UB_LL(DENS)=UB_L(DENS,i,j)
  ! rotate velocity in normal and tangential direction 
  UB_LL(VEL1)=DOT_PRODUCT(UB_L(VELV,i,j),nv(:,i,j))
  UB_LL(VEL2)=DOT_PRODUCT(UB_L(VELV,i,j),t1(:,i,j))
#if PP_dim==3
  UB_LL(VEL3)=DOT_PRODUCT(UB_L(VELV,i,j),t2(:,i,j))
#else
  UB_LL(VEL3)=0.
#endif
  UB_LL(PRES)=UB_L(PRES,i,j)
  UB_LL(SOSP)=UB_L(SOSP,i,j)

  ! right state: U_R
  U_RR(DENS)=U_R(DENS,i,j)
  ! rotate velocity in normal and tangential direction 
  U_RR(VEL1)=DOT_PRODUCT(U_R(VELV,i,j),nv(:,i,j))
  U_RR(VEL2)=DOT_PRODUCT(U_R(VELV,i,j),t1(:,i,j))
#if PP_dim==3
  U_RR(VEL3)=DOT_PRODUCT(U_R(VELV,i,j),t2(:,i,j))
#else
  U_RR(VEL3)=0.
#endif
  U_RR(PRES)=U_R(PRES,i,j)
  ! right baseflow state: UB_R
  UB_RR(DENS)=UB_R(DENS,i,j)
  ! rotate velocity in normal and tangential direction 
  UB_RR(VEL1)=DOT_PRODUCT(UB_R(VELV,i,j),nv(:,i,j))
  UB_RR(VEL2)=DOT_PRODUCT(UB_R(VELV,i,j),t1(:,i,j))
#if PP_dim==3
  UB_RR(VEL3)=DOT_PRODUCT(UB_R(VELV,i,j),t2(:,i,j))
#else
  UB_RR(VEL3)=0.
#endif
  UB_RR(PRES)=UB_R(PRES,i,j)
  UB_RR(SOSP)=UB_R(SOSP,i,j)
  CALL Riemann_loc(U_LL,U_RR,UB_LL,UB_RR,F)

  ! Back Rotate the normal flux into Cartesian direction
  Fout(DENS,i,j)=F(DENS)
  Fout(VELV,i,j)=nv(:,i,j)*F(VEL1) + t1(:,i,j)*F(VEL2) + t2(:,i,j)*F(VEL3)
  Fout(PRES,i,j)=F(PRES)
END DO; END DO
END SUBROUTINE Riemann



!==================================================================================================================================
!> Local Lax-Friedrichs (Rusanov) Riemann solver
!==================================================================================================================================
PURE SUBROUTINE Riemann_LF(U_LL,U_RR,UB_LL,UB_RR,F)
! MODULES
USE MOD_Flux,ONLY:EvalEulerFlux1D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVarBase),INTENT(IN) :: UB_LL,UB_RR
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
REAL    :: F_L(PP_nVar),F_R(PP_nVar)
!==================================================================================================================================
! largest eigenvalue
LambdaMax=MAX(ABS(UB_LL(VEL1)),ABS(UB_RR(VEL1)))+MAX(UB_LL(SOSP),UB_RR(SOSP))
CALL EvalEulerFlux1D(U_LL,UB_LL,F_L)
CALL EvalEulerFlux1D(U_RR,UB_RR,F_R)
F = 0.5*((F_L+F_R) - LambdaMax*(U_RR(:) - U_LL(:)))
END SUBROUTINE Riemann_LF



!==================================================================================================================================
!> Local Flux Vector Splitting Riemann solver
!==================================================================================================================================
PURE SUBROUTINE Riemann_FVS(U_LL,U_RR,UB_LL,UB_RR,F)
! MODULES
USE MOD_Equation_Vars ,ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN)     :: U_LL,U_RR  
REAL,DIMENSION(PP_nVarBase),INTENT(IN) :: UB_LL,UB_RR
REAL,DIMENSION(PP_nVar),INTENT(OUT)    :: F
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                            :: i,j
REAL                               :: UB(PP_nVarBase)
REAL,DIMENSION(PP_nVar)            :: Lambda_plus,Lambda_minus,F_plus,F_minus
REAL                               :: absV,absC1,absC2 
REAL                               :: rho0c0,rho0sc0,srho0sc0,sc02
!==================================================================================================================================
!==================================================================================================================================

! use mean of baseflow states 
UB=0.5*(UB_LL+UB_RR)

! Splitted eigenvalues
absV = ABS(UB(VEL1))
absC1= ABS(UB(VEL1)+UB(SOSP))    ! c1=downstream running wave C1=v+c
absC2= ABS(UB(VEL1)-UB(SOSP))    ! c2=upstream running wave   C1=v-c
Lambda_minus(1:3)=0.5*(UB(VEL1)-absV)
Lambda_minus(4)  =0.5*(UB(VEL1)+UB(SOSP)-absC1)
Lambda_minus(5)  =0.5*(UB(VEL1)-UB(SOSP)-absC2)
Lambda_plus(1:3) =0.5*(UB(VEL1)+absV)
Lambda_plus(4)   =0.5*(UB(VEL1)+UB(SOSP)+absC1)
Lambda_plus(5)   =0.5*(UB(VEL1)-UB(SOSP)+absC2)

rho0c0=UB(DENS)*UB(SOSP)
rho0sc0=UB(DENS)/UB(SOSP)
srho0sc0=1./UB(DENS)/UB(SOSP)
sc02=1./UB(SOSP)/UB(SOSP)
F_plus(1) = Lambda_plus(1)*U_LL(DENS)+0.5*rho0sc0*(Lambda_plus(4)-Lambda_plus(5))*U_LL(VEL1)
F_plus(2) = 0.5*(Lambda_plus(4)+Lambda_plus(5))*U_LL(VEL1) + 0.5*srho0sc0*(Lambda_plus(4)-Lambda_plus(5))*U_LL(PRES)
F_plus(3) = Lambda_plus(2)*U_LL(VEL2) 
F_plus(4) = Lambda_plus(3)*U_LL(VEL3) 
F_plus(5) = 0.5*rho0c0*(Lambda_plus(4)-Lambda_plus(5))*U_LL(VEL1) + 0.5*(Lambda_plus(4)+Lambda_plus(5))*U_LL(PRES)

F_minus(1) = Lambda_minus(1)*U_RR(DENS)+0.5*rho0sc0*(Lambda_minus(4)-Lambda_minus(5))*U_RR(VEL1)
F_minus(2) = 0.5*(Lambda_minus(4)+Lambda_minus(5))*U_RR(VEL1) + 0.5*srho0sc0*(Lambda_minus(4)-Lambda_minus(5))*U_RR(PRES)
F_minus(3) = Lambda_minus(2)*U_RR(VEL2) 
F_minus(4) = Lambda_minus(3)*U_RR(VEL3) 
F_minus(5) = 0.5*rho0c0*(Lambda_minus(4)-Lambda_minus(5))*U_RR(VEL1) + 0.5*(Lambda_minus(4)+Lambda_minus(5))*U_RR(PRES)

F = F_plus+F_minus

END SUBROUTINE Riemann_FVS



!==================================================================================================================================
!> Central flux "Riemann"solver
!==================================================================================================================================
PURE SUBROUTINE Riemann_central(U_LL,U_RR,UB_LL,UB_RR,F)
! MODULES
USE MOD_Flux,ONLY:EvalEulerFlux1D
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVar),INTENT(IN) :: U_LL,U_RR
REAL,DIMENSION(PP_nVarBase),INTENT(IN) :: UB_LL,UB_RR
REAL,DIMENSION(PP_nVar),INTENT(OUT):: F
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL    :: LambdaMax
REAL    :: F_L(PP_nVar),F_R(PP_nVar)
!==================================================================================================================================
CALL EvalEulerFlux1D(U_LL,UB_LL,F_L)
CALL EvalEulerFlux1D(U_RR,UB_RR,F_R)
F = 0.5*(F_L+F_R)
END SUBROUTINE Riemann_central


!==================================================================================================================================
!> Finalize Riemann solver routines
!==================================================================================================================================
SUBROUTINE FinalizeRiemann()
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------!
IMPLICIT NONE
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
END SUBROUTINE FinalizeRiemann

END MODULE MOD_Riemann
