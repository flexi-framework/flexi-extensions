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
!> Routine performing time averaging of variables and the preparation to computing fluctuations
!> For the fluctuations we rely on the fact that \f$ U^{'} U^{'} = \overline{U}^2 - \overline{U^2} \f$
!> The terms computed in this routine are therefore the TimeAvg: \f$ \overline{U} \f$ and 
!> the squared solution denoted by Fluc: \f$ \overline{U^2} \f$ 
!==================================================================================================================================
MODULE MOD_AcSources
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER                        :: nMaxVarAc
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

INTERFACE InitAcSources
  MODULE PROCEDURE InitAcSources
END INTERFACE

INTERFACE FinalizeAcSources
  MODULE PROCEDURE FinalizeAcSources
END INTERFACE

INTERFACE AcSources
  MODULE PROCEDURE AcSources
END INTERFACE

PUBLIC::InitAcSources, FinalizeAcSources, AcSources
!==================================================================================================================================
CONTAINS



!==================================================================================================================================
!> Initializes the time averaging variables and builds map from fluctuation quantities to required time averaged variables
!> (e.g. if VelocityMagnitude fluctuations are to be computed the time averages of the velocity components u,v,w are computed)
!> - only variables specified in the variable list can be averaged
!==================================================================================================================================
SUBROUTINE InitAcSources()
! MODULES
USE MOD_Globals
USE MOD_Preproc
USE MOD_ReadInTools, ONLY: CountOption,GETSTR,GETLOGICAL,GETINT
USE MOD_Mesh_Vars,   ONLY: nElems
USE MOD_AnalyzeEquation_Vars
USE MOD_ACSourceterms_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iVar,iVar2
CHARACTER(LEN=255),ALLOCATABLE :: VarNamesAcIni(:), VarNamesAcList(:)
!==================================================================================================================================
#if !(PARABOLIC)
  CALL CollectiveStop(__STAMP__, &
    'Acoustic source terms require parabolic=T!')
#endif

nVarAc = CountOption('VarNameAc')
IF(nVarAc.EQ.0)THEN
  CALL CollectiveStop(__STAMP__, &
    'No quantities for time averaging have been specified. Please specify quantities or disable time averaging!')
#if FV_ENABLED
ELSE
  CALL CollectiveStop(__STAMP__, &
    'Timeaveraging has not been implemented for FV yet!')
#endif
END IF


! Define variables to be averaged
nMaxVarAc=5
ALLOCATE(VarNamesAcList(nMaxVarAc))
VarNamesAcList(1)  ='PressureTimeDeriv'
VarNamesAcList(2)  ='PressureMatDeriv'
VarNamesAcList(3)  ='LambVecX'
VarNamesAcList(4)  ='LambVecY'
VarNamesAcList(5)  ='LambVecZ'

! Read VarNames from ini file
ALLOCATE(VarNamesAcIni(nVarAc))
DO iVar=1,nVarAc
  VarNamesAcIni(iVar)=GETSTR('VarNameAc')
END DO

nCalcAcSources=GETINT('nCalcAcSources')

! Check which variables have to be calculated and create mappings to global variable index (1:nVarout)
ALLOCATE(CalcAc(nMaxVarAc))
CalcAc=.FALSE.

! check each average from ini file
DO iVar=1,nVarAc
  ! check if avg from ini file is in avg list
  iVar2 = GETMAPBYNAME(VarNamesAcIni(iVar),VarNamesAcIni,nMaxVarAc)
  IF(iVar2.NE.-1)THEN
    CalcAc(iVar2) = .TRUE.
  ELSE
    CALL CollectiveStop(__STAMP__, &
    'Specified varname does not exist: ' // VarNamesAcIni(iVar))
  END IF
END DO


nVarAc=0 ! recount nVarAc
DO iVar=1,nMaxVarAc
  IF(CalcAc(iVar)) nVarAc=nVarAc+1
END DO

! Set indices (iAvg,iFluc) and Varnames
ALLOCATE(VarNamesAcOut(nVarAc))
ALLOCATE(iAc(nMaxVarAc))
! iAc     -> Mapping from VariableList to index in UAc array

VarNamesAcOut(:)=''
nVarAc=0
iAc=0
! Build map for Ac
DO iVar=1,nMaxVarAc
  IF(CalcAc(iVar))THEN
    nVarAc=nVarAc+1
    iAc(iVar)=nVarAc
    VarNamesAcOut(nVarAc) = TRIM(VarNamesAcList(iVar))
  END IF
END DO

! Allocate arrays
ALLOCATE(AcSource(nVarAc,0:PP_N,0:PP_N,0:PP_NZ,nElems))
AcSource=0.

DEALLOCATE(VarNamesAcList,VarNamesAcIni)
END SUBROUTINE InitAcSources


!==================================================================================================================================
!> Return index of string VarName in array VarNameList
!==================================================================================================================================
FUNCTION GETMAPBYNAME(VarName,VarNameList,nVarList)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=*),INTENT(IN)    :: VarName                 !< string to be compared
CHARACTER(LEN=*),INTENT(IN)    :: VarNameList(nVarList)   !< list of strings to be searched
INTEGER,INTENT(IN)             :: nVarList                !< length of list
INTEGER                        :: GETMAPBYNAME            !< index of VarName in VarNameList
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: i
!==================================================================================================================================
GETMAPBYNAME=-1
DO i=1,nVarList
  IF(TRIM(VarName).EQ.TRIM(VarNameList(i)))THEN
    GETMAPBYNAME=i
    RETURN
  END IF
END DO
END FUNCTION



!==================================================================================================================================
!> Compute acoustic source terms
!==================================================================================================================================
SUBROUTINE AcSources(Finalize)
! MODULES
USE MOD_Globals
USE MOD_PreProc
USE MOD_DG_Vars      ,ONLY: U,Ut
#if PARABOLIC
USE MOD_Lifting_Vars ,ONLY: GradUx,GradUy,GradUz
#endif
USE MOD_Mesh_Vars    ,ONLY: MeshFile,nElems
USE MOD_HDF5_Output  ,ONLY: WriteAcSourceTerm
USE MOD_EOS          ,ONLY: ConsToPrim
USE MOD_EOS_Vars     ,ONLY: Kappa,KappaM1
USE MOD_Timedisc_Vars,ONLY: dt,t
USE MOD_AnalyzeEquation_Vars
USE MOD_Mathtools    ,ONLY: CROSS
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
LOGICAL,INTENT(IN)              :: Finalize               !< finalized trapezoidal rule and output file
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                         :: i,j,k,iElem
REAL                            :: tFuture
REAL                            :: tmpVars(nVarAc,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: prim(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: Vorticity(3,0:PP_N,0:PP_N,0:PP_NZ)
REAL                            :: Lamb(3,0:PP_N,0:PP_N,0:PP_NZ)
!----------------------------------------------------------------------------------------------------------------------------------
DO iElem=1,nElems
  IF(CalcAc(1)) &  !'PressureTimeDeriv'
    AcSource(iAc(1),:,:,:,iElem) = KappaM1*(Ut(5,:,:,:,iElem)-1/U(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*Ut(2,:,:,:,iElem)  &
                                           + U(3,:,:,:,iElem)*Ut(3,:,:,:,iElem)  &
                                           + U(4,:,:,:,iElem)*Ut(4,:,:,:,iElem)) &
                                       + 0.5/U(1,:,:,:,iElem)**2*Ut(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*U(2,:,:,:,iElem)   &
                                           + U(3,:,:,:,iElem)*U(3,:,:,:,iElem)   &
                                           + U(4,:,:,:,iElem)*U(4,:,:,:,iElem)))

  IF(CalcAc(2)) &  !'PressureMatDeriv'
    AcSource(iAc(2),:,:,:,iElem) = KappaM1*(Ut(5,:,:,:,iElem)-1/U(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*Ut(2,:,:,:,iElem)  &
                                           + U(3,:,:,:,iElem)*Ut(3,:,:,:,iElem)  &
                                           + U(4,:,:,:,iElem)*Ut(4,:,:,:,iElem)) &
                                       + 0.5/U(1,:,:,:,iElem)**2*Ut(1,:,:,:,iElem)*(  &
                                             U(2,:,:,:,iElem)*U(2,:,:,:,iElem)   &
                                           + U(3,:,:,:,iElem)*U(3,:,:,:,iElem)   &
                                           + U(4,:,:,:,iElem)*U(4,:,:,:,iElem))) &
                                           +(U(2,:,:,:,iElem)*gradUx(5,:,:,:,iElem) &
                                           +U(3,:,:,:,iElem)*gradUy(5,:,:,:,iElem) &
                                           +U(4,:,:,:,iElem)*gradUz(5,:,:,:,iElem) &
                                           )/U(1,:,:,:,iElem)

  IF(CalcAc(3).OR.CalcAc(4).OR.CalcAc(5)) THEN  !'LambVec'
    Vorticity=FillVorticity(gradUx(:,:,:,:,iElem),gradUy(:,:,:,:,iElem),gradUz(:,:,:,:,iElem))
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
      Lamb(:,i,j,k)= -CROSS(Vorticity(:,i,j,k),U(MOMV,i,j,k,iElem)/U(DENS,i,j,k,iElem)) 
    END DO; END DO; END DO
    IF(CalcAc(3)) &
     AcSource(iAc(3),:,:,:,iElem) = Lamb(1,:,:,:)
    IF(CalcAc(4)) &
     AcSource(iAc(4),:,:,:,iElem) = Lamb(2,:,:,:)
    IF(CalcAc(5)) &
     AcSource(iAc(5),:,:,:,iElem) = Lamb(3,:,:,:)
  END IF
END DO ! iElem

! write to file
IF(Finalize)THEN
  tFuture=t+dt
  CALL WriteAcSourceTerm(TRIM(MeshFile),t,tFuture,VarNamesAcOut,AcSource,nVarAc)
END IF

END SUBROUTINE AcSources



FUNCTION FillVorticity(gradUx,gradUy,gradUz) RESULT(Vorticity)
!==================================================================================================================================
! MODULES
USE MOD_Preproc
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
REAL,DIMENSION(PP_nVarPrim,0:PP_N,0:PP_N,0:PP_NZ),INTENT(IN) :: gradUx,gradUy,gradUz
REAL               :: Vorticity(3,0:PP_N,0:PP_N,0:PP_NZ)
!==================================================================================================================================
! VorticityX = dw/dy-dv/dz
  Vorticity(1,:,:,:) = gradUy(4,:,:,:) - gradUz(3,:,:,:)
! VorticityY = du/dz-dw/dx
  Vorticity(2,:,:,:) = gradUz(2,:,:,:) - gradUx(4,:,:,:)
! VorticityZ = dv/dx-du/dy
  Vorticity(3,:,:,:) = gradUx(3,:,:,:) - gradUy(2,:,:,:)
END FUNCTION FillVorticity



!==================================================================================================================================
!> Finalizes the time averaging routines
!==================================================================================================================================
SUBROUTINE FinalizeAcSources()
! MODULES
USE MOD_AnalyzeEquation_Vars
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
SDEALLOCATE(CalcAc)
SDEALLOCATE(iAc)
SDEALLOCATE(AcSource)
SDEALLOCATE(VarNamesAcOut)
END SUBROUTINE FinalizeAcSources

END MODULE MOD_AcSources
