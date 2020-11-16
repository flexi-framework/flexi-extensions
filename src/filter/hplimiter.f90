#if PPLimiter
#if EQNSYSNR == 2
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
MODULE MOD_PPLimiter
!----------------------------------------------------------------------------------------------------------------------------------
! MODULES
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
INTERFACE InitPPLimiter
  MODULE PROCEDURE InitPPLimiter
END INTERFACE

INTERFACE PositivityPreservingLimiter
  MODULE PROCEDURE PositivityPreservingLimiter_Volume
  MODULE PROCEDURE PositivityPreservingLimiter_SidesCons
  MODULE PROCEDURE PositivityPreservingLimiter_SidesPrim
END INTERFACE

INTERFACE PositivityPreservingLimiteriSide
  MODULE PROCEDURE PositivityPreservingLimiteriSide
END INTERFACE

INTERFACE PP_Info
  MODULE PROCEDURE PP_Info
END INTERFACE

PUBLIC:: DefineParametersPPLimiter
PUBLIC:: InitPPLimiter
PUBLIC:: PositivityPreservingLimiter
PUBLIC:: PositivityPreservingLimiteriSide
PUBLIC:: PP_Info
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters needed for filtering
!==================================================================================================================================
SUBROUTINE DefineParametersPPLimiter()
! MODULES
USE MOD_ReadInTools ,ONLY: prms,addStrListEntry
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
CALL prms%SetSection("HP Limiter")
CALL prms%CreateRealOption(            'PPLimiterThreashold',    "Threashold for PP Limiter to modify solution. PPLimiterVar \n"//& 
                                                                 "has to be smaller than this threashold")
CALL prms%CreateRealOption(            'PPLimiterFactor',        "Factor for correction. Is applied to avoid permanent limiting.")                                                           
END SUBROUTINE DefineParametersPPLimiter

!==================================================================================================================================
!> Initialize  information and  operators
!==================================================================================================================================
SUBROUTINE InitPPLimiter()
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_Filter_Vars
USE MOD_ReadInTools        ,ONLY: GETREAL
USE MOD_IO_HDF5            ,ONLY: AddToFieldData
USE MOD_IO_HDF5            ,ONLY: AddToElemData,ElementOut
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

! Read in variables
HPeps = GETREAL('PPLimiterThreashold','1.E-8')
HPfac = GETREAL('PPLimiterFactor','1.00')

! Sanity check
IF (HPfac.GT.1.0) CALL Abort(__STAMP__,'PPLimiterFactor has to be smaller than 1.0!')

! Prepare HP Limiter
ALLOCATE(PP_Elems(nElems))
PP_Elems=0
ALLOCATE(PP_Sides(nElems))
PP_Sides=0
CALL AddToElemData(ElementOut,'HypPresLim',IntArray=PP_Elems)
CALL AddToElemData(ElementOut,'HypPresLimSide',IntArray=PP_Sides)
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
END SUBROUTINE InitPPLimiter

!==================================================================================================================================
!> Hyperbolicity Preserving Limiter, limits polynomial towards admissible cellmean
!==================================================================================================================================
SUBROUTINE PositivityPreservingLimiter_Volume()
! MODULES
USE MOD_PreProc
USE MOD_DG_Vars             ,ONLY: U
USE MOD_Mesh_Vars           ,ONLY: nElems
USE MOD_Filter_Vars         ,ONLY: HPeps,HPfac
USE MOD_Filter_Vars         ,ONLY: PP_Elems,Vol,IntegrationWeight
USE MOD_EOS                 ,ONLY: ConsToPrim,PrimtoCons
#if FV_ENABLED
USE MOD_Filter_Vars         ,ONLY: VolFV,IntegrationWeightFV
USE MOD_FV_Vars             ,ONLY: FV_Elems
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: iElem,i,j,k
REAL                         :: UMean(PP_nVar),rhoMin,pMin
REAL                         :: t,t_loc
!#if FV_RECONSTRUCT
REAL                         :: UPrim(PP_nVarPrim)
!#endif /*FV_RECONSTRUCT*/
!==================================================================================================================================
! mean value
PP_Elems=0.
DO iElem=1,nElems
  UMean = 0.
  rhoMin = HPeps
  pMin = HPeps
  DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
    CALL ConsToPrim(UPrim,U(:,i,j,k,iElem))
    rhoMin=MIN(UPrim(1),rhoMin)
    pMin=MIN(UPrim(5),pMin)
  END DO;END DO;END DO
  IF (rhoMin .GE. HPeps .AND. pMin .GE. HPeps) CYCLE
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
    t = t*HPfac
    DO k=0,PP_NZ;DO j=0,PP_N;DO i=0,PP_N
      U(:,i,j,k,iElem) = t*(U(:,i,j,k,iElem)-UMean) + UMean 
    END DO;END DO;END DO
  END IF
  PP_Elems(iElem)=1
END DO !iElem
END SUBROUTINE PositivityPreservingLimiter_Volume

!==================================================================================================================================
!> Limit if necessary the DG solution at faces between a elements.
!==================================================================================================================================
SUBROUTINE PositivityPreservingLimiter_SidesPrim(UCons_master,UCons_slave,UPrim_master,UPrim_slave)
! MODULES
USE MOD_PreProc
USE MOD_Globals
#if FV_ENABLED
USE MOD_FV_Vars     ,ONLY: FV_Elems_Sum
#endif
USE MOD_Mesh_Vars   ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
USE MOD_Filter_Vars ,ONLY: PP_Sides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: UCons_master(PP_nVar,    0:PP_N,0:PP_NZ,1:nSides) !< Conservative Solution on master side
REAL,INTENT(INOUT)          :: UCons_slave (PP_nVar,    0:PP_N,0:PP_NZ,1:nSides) !< Conservative Solution on slave side
REAL,INTENT(INOUT)          :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Primitive Solution on master side
REAL,INTENT(INOUT)          :: UPrim_slave (PP_nVarPrim,0:PP_N,0:PP_NZ,1:nSides) !< Primitive Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: firstSideID,lastSideID,SideID
REAL        :: UConsTmp_master(PP_nVar,    0:PP_N,0:PP_NZ)              !< 
REAL        :: UConsTmp_slave (PP_nVar,    0:PP_N,0:PP_NZ)              !< 
REAL        :: UPrimTmp_master(PP_nVarPrim,0:PP_N,0:PP_NZ)              !< 
REAL        :: UPrimTmp_slave (PP_nVarPrim,0:PP_N,0:PP_NZ)              !< 
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

PP_Sides = 0
DO SideID=firstSideID,lastSideID
  UConsTmp_slave  = UCons_slave (:,:,:,SideID)
  UConsTmp_master = UCons_master(:,:,:,SideID)
  UPrimTmp_slave  = UPrim_slave (:,:,:,SideID)
  UPrimTmp_master = UPrim_master(:,:,:,SideID)
#if FV_ENABLED
  IF      (FV_Elems_Sum(SideID).EQ.0) THEN
    ! dg
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, UPrimTmp_slave, FVElem=.FALSE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,UPrimTmp_master,FVElem=.FALSE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    ! slave
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, UPrimTmp_slave, FVElem=.TRUE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,UPrimTmp_master,FVElem=.FALSE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.2) THEN
    ! master
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, UPrimTmp_slave, FVElem=.FALSE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,UPrimTmp_master,FVElem=.TRUE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.3) THEN
    ! fv
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave ,UPrimTmp_slave, FVElem=.TRUE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,UPrimTmp_master,FVElem=.TRUE.)
  END IF
#else
  CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave ,UPrimTmp_slave )
  CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,UPrimTmp_master)
#endif
END DO
END SUBROUTINE PositivityPreservingLimiter_SidesPrim

!==================================================================================================================================
!> Limit if necessary the DG solution at faces between a elements.
!==================================================================================================================================
SUBROUTINE PositivityPreservingLimiter_SidesCons(UCons_master,UCons_slave)
! MODULES
USE MOD_PreProc
USE MOD_Globals
#if FV_ENABLED
USE MOD_FV_Vars     ,ONLY: FV_Elems_Sum
#endif
USE MOD_Mesh_Vars   ,ONLY: firstInnerSide,lastMPISide_MINE,nSides
USE MOD_Filter_Vars ,ONLY: PP_Sides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,INTENT(INOUT)          :: UCons_master(PP_nVar,    0:PP_N,0:PP_NZ,1:nSides) !< Conservative Solution on master side
REAL,INTENT(INOUT)          :: UCons_slave (PP_nVar,    0:PP_N,0:PP_NZ,1:nSides) !< Conservative Solution on slave side
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: firstSideID,lastSideID,SideID
REAL        :: UConsTmp_master(PP_nVar,0:PP_N,0:PP_NZ)              !< 
REAL        :: UConsTmp_slave (PP_nVar,0:PP_N,0:PP_NZ)              !< 
!==================================================================================================================================
firstSideID = firstInnerSide
lastSideID  = lastMPISide_MINE

PP_Sides = 0
DO SideID=firstSideID,lastSideID
  UConsTmp_slave  = UCons_slave (:,:,:,SideID)
  UConsTmp_master = UCons_master(:,:,:,SideID)
#if FV_ENABLED
  IF      (FV_Elems_Sum(SideID).EQ.0) THEN
    ! dg
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, FVElem=.FALSE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,FVElem=.FALSE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.1) THEN
    ! slave
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, FVElem=.TRUE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,FVElem=.FALSE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.2) THEN
    ! master
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave, FVElem=.FALSE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,FVElem=.TRUE.)
  ELSE IF (FV_Elems_Sum(SideID).EQ.3) THEN
    ! fv
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave ,FVElem=.TRUE.)
    CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master,FVElem=.TRUE.)
  END IF
#else
  CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_slave )
  CALL PositivityPreservingLimiteriSide(SideID,UConsTmp_master)
#endif
END DO
END SUBROUTINE PositivityPreservingLimiter_SidesCons

!==================================================================================================================================
!> Hyperbolicity Preserving Limiter, limits polynomial towards admissible cellmean
!==================================================================================================================================
SUBROUTINE PositivityPreservingLimiteriSide(SideID,UConsSide,UPrimSide & 
#if FV_ENABLED
  ,FVElem &
#endif
)
! MODULES
USE MOD_PreProc
USE MOD_Filter_Vars         ,ONLY: HPeps,HPfac,PP_Sides
USE MOD_EOS                 ,ONLY: ConsToPrim,PrimtoCons
USE MOD_Mesh_Vars           ,ONLY: SurfElem,SideToElem
USE MOD_Interpolation_Vars  ,ONLY: wGP
#if FV_ENABLED
USE MOD_FV_Vars             ,ONLY: FV_w
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)           :: SideID
REAL,INTENT(INOUT)           :: UConsSide(PP_nVar    ,0:PP_N,0:PP_NZ)
REAL,INTENT(INOUT),OPTIONAL  :: UPrimSide(PP_nVarPrim,0:PP_N,0:PP_NZ)
#if FV_ENABLED
LOGICAL,INTENT(IN),OPTIONAL  :: FVElem
#endif
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                      :: ElemID,i,j
REAL                         :: UMean(PP_nVar),rhoMin,pMin
REAL                         :: t,t_loc
REAL                         :: wloc(0:PP_N)
#if FV_ENABLED
LOGICAL                      :: FVElemloc
#endif
!#if FV_RECONSTRUCT
REAL                         :: UPrimloc(PP_nVarPrim)
REAL                         :: Surf,tmp
!#endif /*FV_RECONSTRUCT*/
!==================================================================================================================================

#if FV_ENABLED
IF (PRESENT(FVElem)) THEN 
  IF (FVElem) wloc(:) = FV_w
ELSE
  wloc(:) = wGP
END IF
#else
  wloc(:) = wGP
#endif
  
! mean value
ElemID = SideToElem(S2E_ELEM_ID,SideID)
IF (ElemID .EQ. -1) RETURN
UMean  = 0.
Surf   = 0.
rhoMin = HPeps
pMin   = HPeps

IF (PRESENT(UPrimSide)) THEN
  rhoMin = MIN(MINVAL(UPrimSide(1,:,:)),rhoMin)
  pMin   = MIN(MINVAL(UPrimSide(5,:,:)),pMin)
ELSE
  DO j=0,PP_NZ;DO i=0,PP_N
    CALL ConsToPrim(UPrimloc,UConsSide(:,i,j))
    rhoMin = MIN(UPrimloc(1),rhoMin)
    pMin   = MIN(UPrimloc(5),pMin)
  END DO;END DO
END IF

IF (rhoMin .GE. HPeps .AND. pMin .GE. HPeps) RETURN
DO j=0,PP_NZ;DO i=0,PP_N
#if PP_dim == 3
  tmp = wloc(i)*wloc(j)*SurfElem(i,j,1,SideID)
#else
  tmp = wloc(i)*SurfElem(i,j,1,SideID)
#endif
  UMean = UMean + UConsSide(:,i,j)*tmp
  Surf  = Surf + tmp 
END DO; END DO
UMean = UMean / Surf
t=1.
DO j=0,PP_NZ;DO i=0,PP_N
  CALL GetTheta(UConsSide(:,i,j),UMean,rhoMin,t_loc)
  t=MIN(t,t_loc)
END DO;END DO
IF(t.LT.1.) THEN
  t = t*HPfac
  DO j=0,PP_NZ;DO i=0,PP_N
    UConsSide(:,i,j) = t*(UConsSide(:,i,j)-UMean) + UMean 
  END DO;END DO
END IF

IF (PRESENT(UPrimSide)) CALL ConsToPrim(PP_N,UPrimSide,UConsSide)
PP_Sides(ElemID)=1
END SUBROUTINE PositivityPreservingLimiteriSide

!==================================================================================================================================
!> Computes thetha, such that theta*U+(1-theta)*cellmean is admissible
!==================================================================================================================================
SUBROUTINE GetTheta(ULoc,UMean,rhoMin,t_out)
! MODULES
USE MOD_PreProc
USE MOD_EOS_Vars      ,ONLY: KappaM1
USE MOD_Filter_Vars   ,ONLY: HPeps
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
REAL                    :: UDiff(PP_nVar)
REAL                    :: p!,p_sGammaM1
!==================================================================================================================================
!p_sGammaM1=ULoc(ENER)-0.5*DOT_PRODUCT(ULoc(MMV2),ULoc(MMV2))/ULoc(DENS)
!IF(ULoc(DENS).GT.HPeps.AND.p_sGammaM1.GT.HPeps) THEN  !TODO: Is HPeps relative or absolute?
!!IF(ULoc(DENS).GT.HPeps) THEN  !TODO: Is HPeps relative or absolute?
  !t_out=1.
  !RETURN  
!END IF


t1 = min((UMean(DENS)-HPeps) / (Umean(DENS)-rhoMin),1.0)
ULoc(DENS)=t1*(ULoc(DENS)-UMean(1))+UMean(1)

p=KappaM1*(ULoc(ENER)-0.5*DOT_PRODUCT(ULoc(MMV2),ULoc(MMV2))/ULoc(DENS))
IF (p .GE. HPeps) THEN
  t_out=1.
  RETURN
END IF
UDiff=ULoc-UMean

a  =                  UDiff(ENER)*UDiff(DENS)  - 0.5*DOT_PRODUCT(UDiff(MMV2),UDiff(MMV2))
b  =      - DOT_PRODUCT(UMean( MMV2),UDiff(MMV2)) + UMean(ENER)*UDiff(DENS) + UMean(DENS)*UDiff(ENER) - (HPeps/KappaM1)*UDiff(DENS)
c  = -0.5*DOT_PRODUCT(UMean( MMV2),UMean( MMV2)) + UMean(ENER)*UMean( DENS) - (HPeps/KappaM1)*UMean(DENS)
IF (a .EQ. 0. .OR. b*b-4.*a*c .LT. 0) THEN 
  t_out = 0.
ELSE
  t_out = -0.5*(b + SQRT(b*b-4.*a*c))/a
END IF
IF(ISNAN(t_out).OR.(t_out.GT.1.).OR.(t_out.LT.0.)) t_out=0.
END SUBROUTINE GetTheta

!==================================================================================================================================
!> Print information on the amount of HP subcells
!==================================================================================================================================
SUBROUTINE PP_Info(iter)
! MODULES
USE MOD_Globals
USE MOD_Mesh_Vars    ,ONLY: nGlobalElems
USE MOD_Analyze_Vars ,ONLY: totalPP_nElems,totalPP_nSides
! IMPLICIT VARIABLE HANDLING
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER(KIND=8),INTENT(IN) :: iter !< number of iterations
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalPP_nElems,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  ! totalPP_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalPP_nElems,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(A,F8.3,A)')' PP amount %: ', totalPP_nElems / REAL(nGlobalElems) / iter*100
totalPP_nElems = 0
#if USE_MPI
IF(MPIRoot)THEN
  CALL MPI_REDUCE(MPI_IN_PLACE,totalPP_nSides,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
  ! totalPP_nElems is counted in PrintStatusLine
ELSE
  CALL MPI_REDUCE(totalPP_nSides,0           ,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_FLEXI,iError)
END IF
#endif
SWRITE(UNIT_stdOut,'(A,F8.3,A)')' PP Sides amount %: ', totalPP_nSides / REAL(nGlobalElems) / iter*100
totalPP_nSides = 0
END SUBROUTINE PP_Info

END MODULE MOD_PPLimiter
#endif
#endif
