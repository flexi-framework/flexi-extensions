#if HPLimiter
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
!==================================================================================================================================
!==================================================================================================================================
MODULE MOD_HPLimiter
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitHPLimiter
  MODULE PROCEDURE InitHPLimiter
END INTERFACE

INTERFACE HyperbolicityPreservingLimiter
  MODULE PROCEDURE HyperbolicityPreservingLimiter
END INTERFACE


PUBLIC:: InitHPLimiter
PUBLIC:: HyperbolicityPreservingLimiter
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Initialize  information and  operators
!==================================================================================================================================
SUBROUTINE InitHPLimiter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_IO_HDF5            ,ONLY: AddToFieldData,FieldOut
USE MOD_Mesh_Vars          ,ONLY: nElems,sJ
USE MOD_Interpolation_Vars ,ONLY: wGP
#if FV_ENABLED
USE MOD_FV_Vars            ,ONLY: FV_w
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k
!==================================================================================================================================
ALLOCATE(t_HPLimiter(1,0:PP_N,0:PP_N,0:PP_NZ,1:nElems))
t_HPLimiter=0.
CALL AddToFieldData(FieldOut,(/1,PP_N+1,PP_N+1,PP_NZ+1/),'HypPresLim',(/'t_HPLimiter'/),RealArray=t_HPLimiter)
IF( .NOT. ALLOCATED(Vol)) THEN 
   ALLOCATE(J_N(0:PP_N,0:PP_N,0:PP_NZ))
   ALLOCATE(Vol(nElems))
   ALLOCATE(Integrationweight(0:PP_N,0:PP_N,0:PP_NZ,nElems))
   Vol = 0.
   DO iElem=1,nElems
     DO k=0,PP_NZ
       DO j=0,PP_N
         DO i=0,PP_N
           J_N(i,j,k)=1./sJ(i,j,k,iElem,0)
#if PP_dim == 3
           IntegrationWeight(i,j,k,iElem) = wGP(i)*wGP(j)*wGP(k)*J_N(i,j,k)
#else
           IntegrationWeight(i,j,k,iElem) = wGP(i)*wGP(j)*J_N(i,j,k)
#endif
           Vol(iElem) = Vol(iElem) + IntegrationWeight(i,j,k,iElem)
         END DO ! i
       END DO ! j
     END DO ! k
   END DO !iElem
   DEALLOCATE(J_N)
END IF
#if FV_ENABLED
   ALLOCATE(J_N(0:PP_N,0:PP_N,0:PP_NZ))
   ALLOCATE(VolFV(nElems))
   ALLOCATE(IntegrationweightFV(0:PP_N,0:PP_N,0:PP_NZ,nElems))
   VolFV = 0.
   DO iElem=1,nElems
     DO k=0,PP_NZ
       DO j=0,PP_N
         DO i=0,PP_N
           J_N(i,j,k)=1./sJ(i,j,k,iElem,1)
#if PP_dim == 3
           IntegrationWeightFV(i,j,k,iElem) = FV_w*FV_w*FV_w*J_N(i,j,k)
#else
           IntegrationWeightFV(i,j,k,iElem) = FV_w*FV_w*J_N(i,j,k)
#endif
           VolFV(iElem) = VolFV(iElem) + IntegrationWeightFV(i,j,k,iElem)
         END DO ! i
       END DO ! j
     END DO ! k
   END DO !iElem
   DEALLOCATE(J_N)
#endif
END SUBROUTINE InitHPLimiter

!==================================================================================================================================
!> Hyperbolicity Preserving Limiter, limits polynomial towards admissible cellmean
!==================================================================================================================================
SUBROUTINE HyperbolicityPreservingLimiter()
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Filter_Vars         ,ONLY: t_HPLimiter,Vol,IntegrationWeight
USE MOD_Interpolation_Vars  ,ONLY: wGP
USE MOD_EOS                 ,ONLY: ConsToPrim,PrimtoCons
#if FV_ENABLED
USE MOD_Filter_Vars         ,ONLY: VolFV,IntegrationWeightFV
USE MOD_FV_Vars             ,ONLY: FV_Elems
!#if FV_RECONSTRUCT
!USE MOD_FV_Vars             ,ONLY: FV_dx_XI_L,FV_dx_ETA_L
!USE MOD_FV_Vars             ,ONLY: FV_dx_XI_R,FV_dx_ETA_R
!USE MOD_FV_Vars             ,ONLY: gradUxi,gradUeta,FV_Elems
#if PP_dim == 3
!USE MOD_FV_Vars             ,ONLY: gradUzeta
!USE MOD_FV_Vars             ,ONLY: FV_dx_ZETA_L
!USE MOD_FV_Vars             ,ONLY: FV_dx_ZETA_R
!#endif
#endif
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k,iVar
REAL                         :: UMean(PP_nVar),rhoMin,pMin
REAL                         :: t,t_loc
REAL                         :: ULoc(PP_nVar)
!#if FV_RECONSTRUCT
REAL                         :: UPrim(PP_nVarPrim)
INTEGER                      :: ii,jj,kk
REAL                         :: UPrim_FV(PP_nVarPrim,0:1,0:1,0:ZDIM(1))
REAL                         :: eps=1.E-8
!#endif /*FV_RECONSTRUCT*/
!==================================================================================================================================
! mean value
DO iElem=1,nElems
  UMean = 0.
  rhoMin = eps
  pMin = eps
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    CALL ConsToPrim(UPrim,U(:,i,j,k,iElem))
    rhoMin=MIN(UPrim(1),rhoMin)
    pMin=MIN(UPrim(5),pMin)
  END DO;END DO;END DO
  IF (rhoMin .GE. eps .AND. pMin .GE. eps) CYCLE
#if FV_ENABLED
  IF (FV_Elems(iElem).EQ.0) THEN
#endif
     DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
       UMean = UMean + U(:,i,j,k,iElem)*IntegrationWeight(i,j,k,iElem)
     END DO; END DO; END DO
     UMean = UMean / Vol(iElem) 
#if FV_ENABLED
  ELSE 
     DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
       UMean = UMean + U(:,i,j,k,iElem)*IntegrationWeightFV(i,j,k,iElem)
     END DO; END DO; END DO
     UMean = UMean / VolFV(iElem) 
  END IF
#endif
  t=1.
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    CALL GetTheta(U(:,i,j,k,iElem),UMean,rhoMin,t_loc)
    t=MIN(t,t_loc)
  END DO;END DO;END DO
  IF(t.LT.1.) THEN
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(:,i,j,k,iElem) = t*(U(:,i,j,k,iElem)-UMean) + UMean 
      !CALL ConsToPrim(UPrim,U(:,i,j,k,iElem))
    END DO;END DO;END DO
  END IF
  t_HPLimiter(:,:,:,:,iElem)=t
END DO !iElem


!TODO: FV reconstruction limiter
!#if FV_RECONSTRUCT
!DO iElem=1,nElems
  !IF(FV_Elems(iElem).EQ.0) CYCLE
  !DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    !t=0.
    !UMean = U(:,i,j,k,iElem)
    !CALL ConsToPrim(UPrim,U(:,i,j,k,iElem))
    !DO iVar=1,PP_nVarPrim
      !UPrim_FV(iVar,:,:,:) = UPrim(iVar)
      !UPrim_FV(iVar,0,:,:) = UPrim_FV(iVar,0,:,:) -   gradUxi(iVar,i,k,j,iElem) *   FV_dx_XI_L(j,k,i,iElem)
      !UPrim_FV(iVar,1,:,:) = UPrim_FV(iVar,1,:,:) +   gradUxi(iVar,i,k,j,iElem) *   FV_dx_XI_R(j,k,i,iElem)
      !UPrim_FV(iVar,:,0,:) = UPrim_FV(iVar,:,0,:) -  gradUeta(iVar,i,k,j,iElem) *  FV_dx_ETA_L(j,k,i,iElem)
      !UPrim_FV(iVar,:,1,:) = UPrim_FV(iVar,:,1,:) +  gradUeta(iVar,i,k,j,iElem) *  FV_dx_ETA_R(j,k,i,iElem)
!#if PP_dim == 3                                                  
      !UPrim_FV(iVar,:,:,0) = UPrim_FV(iVar,:,:,0) - gradUzeta(iVar,i,k,j,iElem) * FV_dx_ZETA_L(j,k,i,iElem)
      !UPrim_FV(iVar,:,:,1) = UPrim_FV(iVar,:,:,1) + gradUzeta(iVar,i,k,j,iElem) * FV_dx_ZETA_R(j,k,i,iElem)
!#endif                                                            
    !END DO
    !DO kk=0,ZDIM(1);DO jj=0,1;DO ii=0,1
      !CALL PrimtoCons(UPrim_FV(:,ii,jj,kk),ULoc)
      !CALL GetTheta(ULoc,UMean,t_loc)
      !t=MAX(t,t_loc)
    !END DO;END DO;END DO
    !IF(t.GT.TINY(1.)) THEN
      !gradUxi  (:,i,k,j,iElem) = (1.-t) *   gradUxi(:,i,k,j,iElem)
      !gradUeta (:,i,k,j,iElem) = (1.-t) *  gradUeta(:,i,k,j,iElem)
!#if PP_dim == 3                                                  
      !gradUzeta(:,i,k,j,iElem) = (1.-t) * gradUzeta(:,i,k,j,iElem)
!#endif                                                            
    !END IF
    !t_HPLimiter(1,i,j,k,iElem)=t
  !END DO;END DO;END DO
!END DO
!#endif /*FV_RECONSTRUCT*/
END SUBROUTINE HyperbolicityPreservingLimiter



!==================================================================================================================================
!> Computes thetha, such that theta*U+(1-theta)*cellmean is admissible
!==================================================================================================================================
SUBROUTINE GetTheta(ULoc,UMean,rhoMin,t_out)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars ,ONLY:KappaM1
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT)      :: ULoc(PP_nVar)
REAL,INTENT(IN)         :: UMean(PP_nVar)
REAL,INTENT(IN)         :: rhoMin
REAL,INTENT(OUT)        :: t_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL                    :: t1,a,b,c
REAL                    :: eps=1.E-8
REAL                    :: UDiff(PP_nVar)
REAL                    :: p!,p_sGammaM1
!==================================================================================================================================
!p_sGammaM1=ULoc(ENER)-0.5*DOT_PRODUCT(ULoc(MMV2),ULoc(MMV2))/ULoc(DENS)
!IF(ULoc(DENS).GT.eps.AND.p_sGammaM1.GT.eps) THEN  !TODO: Is eps relative or absolute?
!!IF(ULoc(DENS).GT.eps) THEN  !TODO: Is eps relative or absolute?
  !t_out=1.
  !RETURN  
!END IF


t1 = min((UMean(DENS)-eps) / (Umean(DENS)-rhoMin),1.0)
ULoc(DENS)=t1*(ULoc(DENS)-UMean(1))+UMean(1)

p=KappaM1*(ULoc(ENER)-0.5*DOT_PRODUCT(ULoc(MMV2),ULoc(MMV2))/ULoc(DENS))
IF (p .GE. eps) THEN
  t_out=1.
  RETURN
END IF
UDiff=ULoc-UMean

a  =                  UDiff(ENER)*UDiff(DENS)  - 0.5*DOT_PRODUCT(UDiff(MMV2),UDiff(MMV2))
b  =      - DOT_PRODUCT(UMean( MMV2),UDiff(MMV2)) + UMean(ENER)*UDiff(DENS) + UMean(DENS)*UDiff(ENER) - (eps/KappaM1)*UDiff(DENS)
c  = -0.5*DOT_PRODUCT(UMean( MMV2),UMean( MMV2)) + UMean(ENER)*UMean( DENS) - (eps/KappaM1)*UMean(DENS)
t_out = -0.5*(b + SQRT(b*b-4.*a*c))/a
IF(ISNAN(t_out).OR.(t_out.GT.1.).OR.(t_out.LT.0.)) t_out=0.
END SUBROUTINE GetTheta

END MODULE MOD_HPLimiter
#endif
