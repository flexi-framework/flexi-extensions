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
!> Contains the GCL volume integral.
!> Computes the volume integral contribution based on mesh velocity and updates Jac_t.
!===================================================================================================================================
MODULE MOD_GCL_VolInt
! MODULES
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
PRIVATE
!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
ABSTRACT INTERFACE
  SUBROUTINE GCL_VolInt_interf(Jac_t)
    USE MOD_PreProc
    USE MOD_Mesh_Vars,ONLY:nElems
    REAL,DIMENSION(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems),INTENT(INOUT) :: Jac_t
  END SUBROUTINE
END INTERFACE

PROCEDURE(GCL_VolInt_interf),POINTER :: GCL_VolInt_pointer   !< pointer defining the GCL volume intergral 
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

PUBLIC::GCL_VolInt_pointer,GCL_VolInt_weakForm
#ifdef SPLIT_DG
PUBLIC::GCL_VolInt_splitForm
#endif /*SPLIT_DG*/
!===================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Computes the volume integral of the weak DG form a la Kopriva.
!> Attention: input Jac_t=0. and is updated with the volume flux derivatives.
!===================================================================================================================================
SUBROUTINE GCL_VolInt_weakForm(Jac_t)
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars        ,ONLY:D_hat_T,nDOFElem
USE MOD_Mesh_Vars      ,ONLY:Metrics_fTilde,Metrics_gTilde,nElems
#if (PP_dim == 3)
USE MOD_Mesh_Vars      ,ONLY:Metrics_hTilde
#endif
USE MOD_MoveMesh_Vars  ,ONLY:Elem_vGP
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
! INPUT VARIABLES
!-----------------------------------------------------------------------------------------------------------------------------------
! OUTPUT VARIABLES
REAL,INTENT(INOUT)                                :: Jac_t(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems)   !< Time derivative of GCL
!-----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_NZ)             :: f,g                  ! volume fluxes at all Gauss points
#if (PP_dim == 3)
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_NZ)             :: h                    ! volume fluxes at all Gauss points
#endif
INTEGER                                           :: i,j,k,iElem          ! Loop indizes
INTEGER                                           :: l                    ! row index for matrix vector product
!===================================================================================================================================
DO iElem=1,nElems
  ! Cut out the local solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components, which consist of the negative mesh velocity
  f(:,:,:) = -Elem_vGP(1,:,:,:,iElem)
  g(:,:,:) = -Elem_vGP(2,:,:,:,iElem)
#if (PP_dim == 3)
  h(:,:,:) = -Elem_vGP(3,:,:,:,iElem)
#endif

  ! Apply transformation
#if (PP_dim == 3)
  CALL GCL_VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem,0),&
                                     Metrics_gTilde(:,:,:,:,iElem,0),&
                                     Metrics_hTilde(:,:,:,:,iElem,0))
#else
  CALL GCL_VolInt_Metrics(nDOFElem,f,g,Metrics_fTilde(:,:,:,:,iElem,0),&
                                     Metrics_gTilde(:,:,:,:,iElem,0))
#endif

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    DO l=0,PP_N
      ! Update the GCL time derivative with the spatial derivatives of the transformed fluxes
      Jac_t(1,i,j,k,iElem) = Jac_t(1,i,j,k,iElem) + D_Hat_T(l,i)*f(l,j,k) + &
#if (PP_dim == 3)
                                                    D_Hat_T(l,k)*h(i,j,l) + &
#endif
                                                    D_Hat_T(l,j)*g(i,l,k)
    END DO ! l
  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE GCL_VolInt_weakForm


#ifdef SPLIT_DG
!==================================================================================================================================
!> Computes the split form of the volume integral of the weak DG form.
!> Attention: input Jac_t=0. and is updated with the volume flux derivatives.
!==================================================================================================================================
SUBROUTINE GCL_VolInt_splitForm(Jac_t)
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: DVolSurf,nDOFElem
USE MOD_Mesh_Vars    ,ONLY: Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
USE MOD_MoveMesh_Vars,ONLY: Elem_vGP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)   :: Jac_t(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems) !< Time derivative of Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                               :: i,j,k,l,iElem
REAL                                  :: Flux   !< temp variable for split flux
REAL,DIMENSION(0:PP_N,0:PP_N,0:PP_NZ) :: f,g,h  !< GCL fluxes at GP
!==================================================================================================================================
! Diffusive part
DO iElem=1,nElems
  ! Cut out the local solution for a grid cell iElem and all Gauss points from the global field
  ! Compute for all Gauss point values the Cartesian flux components, which consist of the negative mesh velocity
  f = -Elem_vGP(1,:,:,:,iElem)
  g = -Elem_vGP(2,:,:,:,iElem)
#if (PP_dim == 3)
  h = -Elem_vGP(3,:,:,:,iElem)
#endif
#if (PP_dim == 3)
  CALL GCL_VolInt_Metrics(nDOFElem,f,g,h,Metrics_fTilde(:,:,:,:,iElem,0),&
                                         Metrics_gTilde(:,:,:,:,iElem,0),&
                                         Metrics_hTilde(:,:,:,:,iElem,0))
#else
  CALL GCL_VolInt_Metrics(nDOFElem,f,g,Metrics_fTilde(:,:,:,:,iElem,0),&
                                       Metrics_gTilde(:,:,:,:,iElem,0))
#endif

  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N

    DO l=i+1,PP_N
       ! compute split flux in x-direction
       Flux = -0.5 * ((Elem_vGP(1,i,j,k,iElem) + Elem_vGP(1,l,j,k,iElem)) * (Metrics_fTilde(1,i,j,k,iElem,0) + Metrics_fTilde(1,l,j,k,iElem,0)) + &
                      (Elem_vGP(2,i,j,k,iElem) + Elem_vGP(2,l,j,k,iElem)) * (Metrics_fTilde(2,i,j,k,iElem,0) + Metrics_fTilde(2,l,j,k,iElem,0)) + &
                      (Elem_vGP(3,i,j,k,iElem) + Elem_vGP(3,l,j,k,iElem)) * (Metrics_fTilde(3,i,j,k,iElem,0) + Metrics_fTilde(3,l,j,k,iElem,0)))
       Jac_t(:,i,j,k,iElem) = Jac_t(:,i,j,k,iElem) + DVolSurf(l,i) * Flux
       Jac_t(:,l,j,k,iElem) = Jac_t(:,l,j,k,iElem) + DVolSurf(i,l) * Flux
    END DO ! m

    DO l=j+1,PP_N
       ! compute split flux in y-direction
       Flux = -0.5 * ((Elem_vGP(1,i,j,k,iElem) + Elem_vGP(1,i,l,k,iElem)) * (Metrics_gTilde(1,i,j,k,iElem,0) + Metrics_gTilde(1,i,l,k,iElem,0)) + &
                      (Elem_vGP(2,i,j,k,iElem) + Elem_vGP(2,i,l,k,iElem)) * (Metrics_gTilde(2,i,j,k,iElem,0) + Metrics_gTilde(2,i,l,k,iElem,0)) + &
                      (Elem_vGP(3,i,j,k,iElem) + Elem_vGP(3,i,l,k,iElem)) * (Metrics_gTilde(3,i,j,k,iElem,0) + Metrics_gTilde(3,i,l,k,iElem,0)))
       Jac_t(:,i,j,k,iElem) = Jac_t(:,i,j,k,iElem) + DVolSurf(l,j) * Flux
       Jac_t(:,i,l,k,iElem) = Jac_t(:,i,l,k,iElem) + DVolSurf(j,l) * Flux
    END DO ! m

#if PP_dim==3
    DO l=k+1,PP_N
       ! compute split flux in z-direction
       Flux = -0.5 * ((Elem_vGP(1,i,j,k,iElem) + Elem_vGP(1,i,j,l,iElem)) * (Metrics_hTilde(1,i,j,k,iElem,0) + Metrics_hTilde(1,i,j,l,iElem,0)) + &
                      (Elem_vGP(2,i,j,k,iElem) + Elem_vGP(2,i,j,l,iElem)) * (Metrics_hTilde(2,i,j,k,iElem,0) + Metrics_hTilde(2,i,j,l,iElem,0)) + &
                      (Elem_vGP(3,i,j,k,iElem) + Elem_vGP(3,i,j,l,iElem)) * (Metrics_hTilde(3,i,j,k,iElem,0) + Metrics_hTilde(3,i,j,l,iElem,0)))
       Jac_t(:,i,j,k,iElem) = Jac_t(:,i,j,k,iElem) + DVolSurf(l,k) * Flux
       Jac_t(:,i,j,l,iElem) = Jac_t(:,i,j,l,iElem) + DVolSurf(k,l) * Flux
    END DO ! l
#endif /*PP_dim==3*/


    ! consistency: if both points for flux evaluation are the same the standart
    ! euler fluxes are retained
    Jac_t(:,i,j,k,iElem) = Jac_t(:,i,j,k,iElem) + DVolSurf(i,i)*f(i,j,k) + &
#if PP_dim==3
                                                  DVolSurf(k,k)*h(i,j,k) + &
#endif /*PP_dim==3*/
                                                  DVolSurf(j,j)*g(i,j,k)

  END DO; END DO; END DO !i,j,k
END DO ! iElem
END SUBROUTINE GCL_VolInt_splitForm
#endif /*SPLIT_DG*/

!==================================================================================================================================
!> Compute the tranformed states for all conservative variables using the metric terms
!==================================================================================================================================
#if (PP_dim == 3)
SUBROUTINE GCL_VolInt_Metrics(nDOFs,f,g,h,Mf,Mg,Mh)
#else
SUBROUTINE GCL_VolInt_Metrics(nDOFs,f,g,  Mf,Mg)
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
USE MOD_PreProc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)                          :: nDOFs    !< Number of DOFs per element
                                                        !> Metric terms
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mf,Mg
#if (PP_dim == 3)
REAL,DIMENSION(3,nDOFs),INTENT(IN)          :: Mh
#endif
                                                        !> Volume fluxes at GP to be transformed from physical to reference space
REAL,DIMENSION(1,nDOFs),INTENT(INOUT)       :: f,g
#if (PP_dim == 3)
REAL,DIMENSION(1,nDOFs),INTENT(INOUT)       :: h
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                                     :: i
REAL,DIMENSION(1)                           :: fTilde,gTilde !< Auxiliary variables needed to store the fluxes at one GP
#if (PP_dim == 3)
REAL,DIMENSION(1)                           :: hTilde        !< Auxiliary variables needed to store the fluxes at one GP
#endif
!==================================================================================================================================
DO i=1,nDOFs
  fTilde=f(:,i)
  gTilde=g(:,i)
#if (PP_dim == 3)
  hTilde=h(:,i)
  ! Compute the transformed fluxes with the metric terms
  ! Attention 1: we store the transformed fluxes in f,g,h again
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i) + &
           hTilde*Mf(3,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i) + &
           hTilde*Mg(3,i)
  h(:,i) = fTilde*Mh(1,i) + &
           gTilde*Mh(2,i) + &
           hTilde*Mh(3,i)
#else
  f(:,i) = fTilde*Mf(1,i) + &
           gTilde*Mf(2,i)
  g(:,i) = fTilde*Mg(1,i) + &
           gTilde*Mg(2,i)
#endif
END DO ! i
END SUBROUTINE GCL_VolInt_Metrics

END MODULE MOD_GCL_VolInt
