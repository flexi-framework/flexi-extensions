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
!> This module contains routines necessary for the Equation of State. For the linear scalar advection diffusion equation,
!> no EOS is needed. This module provides some dummy routines due to compatibility issues with the Navier-Stokes equation system.
!==================================================================================================================================
MODULE MOD_Baseflow
! MODULES
IMPLICIT NONE
PRIVATE

INTEGER,PARAMETER :: BASEFLOWTYPE_CONST      = 0
INTEGER,PARAMETER :: BASEFLOWTYPE_ANALYTIC   = 1
INTEGER,PARAMETER :: BASEFLOWTYPE_FILE       = 2
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! Private Part ---------------------------------------------------------------------------------------------------------------------
! Public Part ----------------------------------------------------------------------------------------------------------------------

INTERFACE InitBaseflow
  MODULE PROCEDURE InitBaseflow
END INTERFACE

INTERFACE CalcSourceBF
  MODULE PROCEDURE CalcSourceBF
END INTERFACE

PUBLIC::DefineParametersBaseflow, InitBaseflow, CalcSourceBF, FinalizeBaseflow
!==================================================================================================================================

CONTAINS

!==================================================================================================================================
!> Define parameters 
!==================================================================================================================================
SUBROUTINE DefineParametersBaseflow()
! MODULES
USE MOD_ReadInTools ,ONLY: prms
USE MOD_ReadInTools ,ONLY: addStrListEntry
IMPLICIT NONE
!==================================================================================================================================
CALL prms%SetSection("AcousticBaseFlow")
CALL prms%CreateIntFromStringOption('BaseflowType', "Type of baseflow for acoustic equation. Constant (0), analytic (1), file (2)",&
                                                    'constant')
CALL addStrListEntry('baseFlowType','constant',      BASEFLOWTYPE_CONST)
CALL addStrListEntry('baseFlowType','analytic',      BASEFLOWTYPE_ANALYTIC)
CALL addStrListEntry('baseFlowType','file',          BASEFLOWTYPE_FILE)
CALL prms%CreateRealArrayOption('BaseState',     "Background State in primitive variables (density, velx, vely, velz, pressure).")
CALL prms%CreateIntOption('BaseflowFunc',    "Analytical baseflow function. 1: potential flow cylinder, 2: shear-flow.")
CALL prms%CreateStringOption( 'BaseFlowFile',"FLEXI solution (e.g. TimeAvg) file from which baseflow is read.")
CALL prms%CreateLogicalOption('Baseflow_visu',   "Output a visualization of the baseflow state in vtu format.", 'F')
CALL prms%CreateRealOption(   'shearlayerPeakVelocity', "Peak velocity of the shear-layer Analytical baseflow 2")
CALL prms%CreateRealOption(   'shearlayerThickness',    "thickness of the shear-layer Analytical baseflow 2")
CALL prms%CreateRealOption(   'mixinglayerThickness',    "thickness of the mixing-layer Analytical baseflow 3")
CALL prms%CreateRealOption(   'u1',    "... baseflow 3")
CALL prms%CreateRealOption(   'u2',    "... baseflow 3")
END SUBROUTINE DefineParametersBaseflow


!==================================================================================================================================
!> Initialize the base flow and its gradients for the linearized Euler equations
!==================================================================================================================================
SUBROUTINE InitBaseflow()
! MODULES
USE MOD_ReadInTools,   ONLY: GETINTFROMSTR,GETREALARRAY,GETSTR,GETINT,GETLOGICAL,GETREAL
USE MOD_Preproc
USE MOD_Globals
USE MOD_Equation_Vars, ONLY:Kappa,StrVarNamesBase
USE MOD_Baseflow_Vars
USE MOD_Mesh_Vars,     ONLY: nElems,nSides,Elem_xGP
USE MOD_Output_Vars,   ONLY: NVisu,Vdm_GaussN_NVisu,ProjectName
USE MOD_ChangeBasis,   ONLY: ChangeBasis3D
USE MOD_VTK           ,ONLY:WriteDataToVTK
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(LEN=300) :: baseFlowFile,FileString
INTEGER            :: baseflowFunc
LOGICAL            :: Baseflow_visu
INTEGER            :: iElem
REAL,ALLOCATABLE   :: UBase_NVisu(:,:,:,:,:),Coords_NVisu(:,:,:,:,:)
!==================================================================================================================================
baseFlowType=GETINTFROMSTR('BaseFlowType') ! constant (0), analytic (1), file (2)
! Read in base state
BaseState(1:5) = GETREALARRAY('BaseState',5)
BaseState(6)=SQRT(Kappa*BaseState(5)/BaseState(1))
varMeanFlow=.FALSE.
SELECT CASE(BaseFlowType)
CASE(BASEFLOWTYPE_CONST)
CASE(BASEFLOWTYPE_ANALYTIC) ! analytic, specify which one...
  baseflowFunc=GETINT('baseflowFunc')
  SELECT CASE(baseflowFunc)
  CASE(2)
    uShear=GETREAL('shearlayerPeakVelocity')
    dShear=GETREAL('shearlayerThickness')
  CASE(3)
    dShear=GETREAL('mixinglayerThickness')
    u1=GETREAL('u1')
    u2=GETREAL('u2')
  END SELECT

varMeanFlow=.TRUE.
CASE(BASEFLOWTYPE_FILE)     ! file, specify base flow file.
  baseflowFile=GETSTR('baseflowFile')
varMeanFlow=.TRUE.
END SELECT
baseFlow_visu=GETLOGICAL('baseFlow_visu','F')

! Allocate variable base flow
ALLOCATE(UBase(1:PP_nVarBase,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUbx(1:PP_nVarBase,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUby(1:PP_nVarBase,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(gradUbz(1:PP_nVarBase,0:PP_N,0:PP_N,0:PP_NZ,nElems))
ALLOCATE(UBase_Master(PP_nVarBase,0:PP_N,0:PP_NZ,1:nSides))
ALLOCATE(UBase_Slave(PP_nVarBase,0:PP_N,0:PP_NZ,1:nSides))

! Initialize Base flow
SELECT CASE(BaseFlowType)
CASE(BASEFLOWTYPE_CONST)
UBase(1,:,:,:,:)=BaseState(1)
UBase(2,:,:,:,:)=BaseState(2)
UBase(3,:,:,:,:)=BaseState(3)
UBase(4,:,:,:,:)=BaseState(4)
UBase(5,:,:,:,:)=BaseState(5)
UBase(6,:,:,:,:)=BaseState(6)
! TODO : rotate baseflow surface states to local side system to avoid extra runtime cost
UBase_Slave(1,:,:,:)=BaseState(1)
UBase_Slave(2,:,:,:)=BaseState(2)
UBase_Slave(3,:,:,:)=BaseState(3)
UBase_Slave(4,:,:,:)=BaseState(4)
UBase_Slave(5,:,:,:)=BaseState(5)
UBase_Slave(6,:,:,:)=BaseState(6)
UBase_Slave(1,:,:,:)=BaseState(1)
UBase_Master(1,:,:,:)=BaseState(1)
UBase_Master(2,:,:,:)=BaseState(2)
UBase_Master(3,:,:,:)=BaseState(3)
UBase_Master(4,:,:,:)=BaseState(4)
UBase_Master(5,:,:,:)=BaseState(5)
UBase_Master(6,:,:,:)=BaseState(6)
gradUbx=0.
gradUby=0.
gradUbz=0.

CASE(BASEFLOWTYPE_ANALYTIC)
  CALL FillAnalyticBaseflow(baseflowFunc,UBase)

CASE(BASEFLOWTYPE_FILE)
  CALL ReadBaseflow(baseflowFile)
END SELECT

! Calculate gradients and surface distributions of the baseflow
IF(varMeanFlow) CALL CalcBaseGradandSurf()

! Baseflow visualization
IF(Baseflow_visu) THEN
  FileString=TRIM(INTSTAMP(TRIM(ProjectName),myRank))//'_Baseflow.vtu'
  ALLOCATE(Coords_NVisu(1:3, 0:NVisu,0:NVisu,0:NVisu,nElems))
  ALLOCATE(UBase_NVisu(PP_nVarBase,0:NVisu,0:NVisu,0:NVisu,nElems))
  ! Create coordinates of visualization points
  DO iElem=1,nElems
    CALL ChangeBasis3D(3,PP_N,NVisu,Vdm_GaussN_NVisu,Elem_xGP(1:3,:,:,:,iElem),Coords_NVisu(1:3,:,:,:,iElem))
  END DO
  ! Interpolate solution onto visu grid
  UBase_NVisu=0.
  DO iElem=1,nElems
    CALL ChangeBasis3D(PP_nVarBase,PP_N,NVisu,Vdm_GaussN_NVisu,UBase(:,:,:,:,iElem),UBase_NVisu(:,:,:,:,iElem))
  END DO !SpongeElem=1,nSpongeElems
  CALL WriteDataToVTK(PP_nVarBase,NVisu,nElems,StrVarNamesBase,Coords_NVisu,UBase_NVisu,TRIM(FileString),dim=PP_dim)
  DEALLOCATE(Coords_NVisu)
  DEALLOCATE(UBase_NVisu)
END IF !Baseflow_visu

END SUBROUTINE InitBaseflow


!==================================================================================================================================
!> Compute source terms of LEE arising from mean flow gradients. Attention: this routine is also used to nullify the source array!
!==================================================================================================================================
SUBROUTINE CalcSourceBF(source)
! MODULES
USE MOD_Equation_Vars, ONLY:Kappa,KappaM1
USE MOD_Baseflow_Vars, ONLY:UBase,gradUbx,gradUby,gradUbz
USE MOD_Mesh_Vars,     ONLY:nElems
USE MOD_PreProc
USE MOD_DG_Vars,       ONLY:U
USE MOD_Equation_Vars, ONLY: IniExactFunc,IniRefState
USE MOD_Mesh_Vars,     ONLY:Elem_xGP
USE MOD_Baseflow_Vars, ONLY:uShear,dShear
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(OUT)  :: source(PP_nVar,0:PP_N,0:PP_N,0:PP_NZ,nElems)  !< source array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                  :: iElem,i,j,k,l
REAL                     :: divV,srho
REAL,DIMENSION(PP_dim)   :: gradub,gradvb,gradwb,gradsrho,gradp
!REAL,DIMENSION(PP_nVar,PP_nVar)  :: D,Dx
!==================================================================================================================================
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
#if PP_dim==3
    divV = gradUbx(2,i,j,k,iElem)+gradUby(3,i,j,k,iElem)+gradUbz(4,i,j,k,iElem)
    srho = 1./UBase(1,i,j,k,iElem)

    gradub  =(/gradUbx(2,i,j,k,iElem),gradUby(2,i,j,k,iElem),gradUbz(2,i,j,k,iElem)/)
    gradvb  =(/gradUbx(3,i,j,k,iElem),gradUby(3,i,j,k,iElem),gradUbz(3,i,j,k,iElem)/)
    gradwb  =(/gradUbx(4,i,j,k,iElem),gradUby(4,i,j,k,iElem),gradUbz(4,i,j,k,iElem)/)
    gradsrho=-srho**2*(/gradUbx(1,i,j,k,iElem),gradUby(1,i,j,k,iElem),gradUbz(1,i,j,k,iElem)/) ! (1/rho_0)_x=-1/rho_0^2*(rho_0)_x
    gradp   =(/gradUbx(5,i,j,k,iElem),gradUby(5,i,j,k,iElem),gradUbz(5,i,j,k,iElem)/)


    source(DENS,i,j,k,iElem)=0.
    source(VEL1,i,j,k,iElem)=-srho*DOT_PRODUCT(UBase(VELV,i,j,k,iELem),gradub)*U(DENS,i,j,k,iElem) &
                                                              + (divV-gradub(1))*U(VEL1,i,j,k,iElem) &
                                                              -       gradub(2) *U(VEL2,i,j,k,iElem) &
                                                              -       gradub(3) *U(VEL3,i,j,k,iElem) &
                                                              +     gradsrho(1) *U(PRES,i,j,k,iElem)

    source(VEL2,i,j,k,iElem)=-srho*DOT_PRODUCT(UBase(VELV,i,j,k,iELem),gradvb)*U(DENS,i,j,k,iElem) &
                                                              -       gradvb(1) *U(VEL1,i,j,k,iElem) &
                                                              + (divV-gradvb(2))*U(VEL2,i,j,k,iElem) &
                                                              -       gradvb(3) *U(VEL3,i,j,k,iElem) &
                                                              +     gradsrho(2) *U(PRES,i,j,k,iElem)

    source(VEL3,i,j,k,iElem)=-srho*DOT_PRODUCT(UBase(VELV,i,j,k,iELem),gradwb)*U(DENS,i,j,k,iElem) &
                                                              -       gradwb(1) *U(VEL1,i,j,k,iElem) &
                                                              -       gradwb(2) *U(VEL2,i,j,k,iElem) &
                                                              + (divV-gradwb(3))*U(VEL3,i,j,k,iElem) &
                                                              +     gradsrho(3) *U(PRES,i,j,k,iElem)

    source(PRES,i,j,k,iElem)=                                 KappaM1*gradp(1) *U(VEL1,i,j,k,iElem) &
                                                              +KappaM1*gradp(2) *U(VEL2,i,j,k,iElem) &
                                                              +KappaM1*gradp(3) *U(VEL3,i,j,k,iElem) &
                                                              -KappaM1*divV     *U(PRES,i,j,k,iElem)
#else
! 2D variant
    divV = gradUbx(2,i,j,k,iElem)+gradUby(3,i,j,k,iElem)
    srho = 1./UBase(1,i,j,k,iElem)

    gradub  =(/gradUbx(2,i,j,k,iElem),gradUby(2,i,j,k,iElem)/)
    gradvb  =(/gradUbx(3,i,j,k,iElem),gradUby(3,i,j,k,iElem)/)
    gradsrho=-srho**2*(/gradUbx(1,i,j,k,iElem),gradUby(1,i,j,k,iElem)/) ! (1/rho_0)_x=-1/rho_0^2*(rho_0)_x
    gradp   =(/gradUbx(5,i,j,k,iElem),gradUby(5,i,j,k,iElem)/)


    source(DENS,i,j,k,iElem)=0.
    source(VEL1,i,j,k,iElem)=-srho*DOT_PRODUCT(UBase(2:3,i,j,k,iELem),gradub)*U(DENS,i,j,k,iElem) &

                                                              + (divV-gradub(1))*U(VEL1,i,j,k,iElem) &
                                                              -       gradub(2) *U(VEL2,i,j,k,iElem) &
                                                              +     gradsrho(1) *U(PRES,i,j,k,iElem)

    source(VEL2,i,j,k,iElem)=-srho*DOT_PRODUCT(UBase(2:3,i,j,k,iELem),gradvb)*U(DENS,i,j,k,iElem) &
                                                              -       gradvb(1) *U(VEL1,i,j,k,iElem) &
                                                              + (divV-gradvb(2))*U(VEL2,i,j,k,iElem) &
                                                              +     gradsrho(2) *U(PRES,i,j,k,iElem)

    source(VEL3,i,j,k,iElem)=0.
    source(PRES,i,j,k,iElem)=                                 KappaM1*gradp(1) *U(VEL1,i,j,k,iElem) &
                                                              +KappaM1*gradp(2) *U(VEL2,i,j,k,iElem) &
                                                              -KappaM1*divV     *U(PRES,i,j,k,iElem)

#endif
  END DO; END DO; END DO ! i,j,k
END DO ! iElem
END SUBROUTINE CalcSourceBF


!==================================================================================================================================
!> Compute source terms of LEE arising from mean flow gradients 
!==================================================================================================================================
SUBROUTINE FillAnalyticBaseflow(baseflowFunc,UBase)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Baseflow_Vars, ONLY:BaseState,uShear,dShear,u1,u2
USE MOD_Equation_Vars, ONLY:Kappa
USE MOD_Mesh_Vars,     ONLY:nElems,Elem_xGP
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)  :: baseflowFunc
REAL,INTENT(INOUT)  :: UBase(PP_nVarBase,0:PP_N,0:PP_N,0:PP_NZ,nElems)  !< source array
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iElem,i,j,k,l
REAL                :: r,uinf,pinf,x(3)
REAL                :: dShearTmp
!==================================================================================================================================
SELECT CASE(baseflowFunc)
CASE(1) ! cylinder
  r=1.
  uinf=BaseState(VEL1)
  pinf=BaseState(PRES)
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    x=Elem_xGP(:,i,j,k,iElem)
    UBase(1,i,j,k,iElem)=1.
    UBase(2,i,j,k,iElem)=uinf*(1.-(r**2*((x(1)**2-x(2)**2)/((x(1)**2+x(2)**2)**2))))
    UBase(3,i,j,k,iElem)=-2.*uinf*r**2*((x(1)*x(2))/((x(1)**2+x(2)**2)**2))
    UBase(4,i,j,k,iElem)=0.
    UBase(5,i,j,k,iElem)=pinf+UBase(1,i,j,k,iElem)/2*(uinf**2-UBase(2,i,j,k,iElem)**2)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE(2) !shearlayer
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    x=Elem_xGP(:,i,j,k,iElem)
    UBase(1,i,j,k,iElem)=BaseState(DENS)
    UBase(2,i,j,k,iElem)=uShear*tanh(2*x(2)/dShear)
    UBase(3,i,j,k,iElem)=0
    UBase(4,i,j,k,iElem)=0
    UBase(5,i,j,k,iElem)=BaseState(PRES)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem
CASE(3) !mixinglayer
  DO iElem=1,nElems
    DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    x=Elem_xGP(:,i,j,k,iElem)
    dShearTmp=dShear*(3/2+1/2*tanh((x(1)-70)/200))
    UBase(1,i,j,k,iElem)=BaseState(DENS)
    UBase(2,i,j,k,iElem)=(u1+u2)/2+((u2-u1)/2)*tanh(2*x(2)/dShearTmp)
    UBase(3,i,j,k,iElem)=0
    UBase(4,i,j,k,iElem)=0
    UBase(5,i,j,k,iElem)=BaseState(PRES)
    END DO; END DO; END DO ! i,j,k
  END DO ! iElem

CASE DEFAULT
  CALL CollectiveStop(__STAMP__,&
    "Baseflowfunc does not exist!")
END SELECT
UBase(6,:,:,:,:)=SQRT(Kappa*UBase(5,:,:,:,:)/UBase(1,:,:,:,:))
END SUBROUTINE FillAnalyticBaseflow


!==================================================================================================================================
!> \brief Read baseflow from HDF5 file
!> 
!> This routine reads the base flow  from a .h5 file. It is checked if the base flow is using the same polynomial degree and
!> node type as the current solution.
!> If not, a interpolation is performed first.
!==================================================================================================================================
SUBROUTINE ReadBaseFlow(FileName)
! MODULES
USE MOD_PreProc
USE MOD_Globals
USE MOD_Baseflow_Vars,     ONLY: UBase
USE MOD_Mesh_Vars,         ONLY: offsetElem,nGlobalElems,nElems
USE MOD_HDF5_input,        ONLY: OpenDataFile,CloseDataFile,ReadArray,GetDataProps,ReadAttribute
USE MOD_IO_HDF5,           ONLY: File_ID
USE MOD_ChangeBasis,       ONLY: ChangeBasis3D
USE MOD_Interpolation,     ONLY: GetVandermonde
USE MOD_Interpolation_Vars,ONLY: NodeType
USE MOD_Equation_Vars,     ONLY: Kappa
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
CHARACTER(LEN=300),INTENT(IN)  :: FileName                 !< HDF5 filename
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                        :: iElem,iVar1,iVar2
INTEGER                        :: N_HDF5,nVar_HDF5,nElems_HDF5
CHARACTER(LEN=255)             :: NodeType_HDF5
CHARACTER(LEN=255),ALLOCATABLE :: VarNames_HDF5(:)
REAL,ALLOCATABLE               :: UTmp(:,:,:,:,:),UTmp2(:,:,:,:,:),Vdm_NHDF5_N(:,:)
INTEGER,PARAMETER              :: nVarCons=5
INTEGER,PARAMETER              :: nVarPrim=5
INTEGER                        :: mapPrim(5),mapCons(5)
INTEGER                        :: N_HDF5Z
LOGICAL                        :: primVars
REAL                           :: KappaM1
CHARACTER(LEN=255),DIMENSION(5),PARAMETER :: StrVarNamesCons =&
  (/ CHARACTER(LEN=255) :: 'Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'/) !< conservative variable names
CHARACTER(LEN=255),DIMENSION(6),PARAMETER :: StrVarNamesPrim=&
  (/ CHARACTER(LEN=255) :: 'Density','VelocityX','VelocityY','VelocityZ','Pressure','Temperature'/) !< primitive variable names
!==================================================================================================================================
SWRITE(UNIT_StdOut,'(A,A)')'  Read Base Flow from file "',FileName
CALL OpenDataFile(FileName,create=.FALSE.,single=.FALSE.,readOnly=.TRUE.)
CALL GetDataProps(nVar_HDF5,N_HDF5,nElems_HDF5,NodeType_HDF5,'Mean')

!IF((nElems_HDF5.NE.nGlobalElems).OR.(nVar_HDF5.NE.PP_nVar))THEN
IF((nElems_HDF5.NE.nGlobalElems).OR.(nVar_HDF5.LT.PP_nVar))THEN
  CALL abort(__STAMP__,&
             'Baseflow file does not match solution. Elements,nVar',nElems_HDF5,REAL(nVar_HDF5))
ENDIF

! check the variable names
ALLOCATE(VarNames_HDF5(nVar_HDF5))
!CALL ReadAttribute(File_ID,'VarNames',nVar_HDF5,StrArray=VarNames_HDF5)
CALL ReadAttribute(File_ID,'VarNames_Mean',nVar_HDF5,StrArray=VarNames_HDF5)
! case 1: primitive. 
mapPrim=-999
primVars=.TRUE.
DO iVar1=1,nVarPrim
  DO iVar2=1,nVar_HDF5
    IF(TRIM(VarNames_HDF5(iVar2)).EQ.TRIM(StrVarNamesPrim(iVar1)))THEN
      mapPrim(iVar1)=iVar2
    END IF
  END DO
END DO
! case 2: conservative state => convert to primitive
IF(.NOT.ANY(mapPrim.EQ.-999))THEN
  SWRITE(UNIT_StdOut,'(A)')' Found complete primitive state vector in baseflow file'
ELSE
  primVars=.FALSE.
  mapCons=-999
  DO iVar1=1,nVarCons
    DO iVar2=1,nVar_HDF5
          IF(TRIM(VarNames_HDF5(iVar2)).EQ.TRIM(StrVarNamesCons(iVar1)))THEN
        mapCons(iVar1)=iVar2
      END IF
    END DO
  END DO
  IF(ANY(mapCons.EQ.-999))THEN
    CALL CollectiveStop(__STAMP__,&
      "Baseflowfile does neither contain a complete set of conservative nor primitive variables!")
  ELSE
    SWRITE(UNIT_StdOut,'(A)')' Found complete conservative state vector in baseflow file. Converting to primitive..'
  END IF
END IF

! Read in state
ALLOCATE(UTmp(nVar_HDF5,0:PP_N,0:PP_N,0:PP_NZ,nElems))
IF((N_HDF5.EQ.PP_N).AND.(TRIM(NodeType_HDF5).EQ.TRIM(NodeType)))THEN
  ! No interpolation needed, read solution directly from file
  CALL ReadArray('Mean',5,(/nVar_HDF5,PP_N+1,PP_N+1,PP_NZ+1,nElems/),OffsetElem,5,RealArray=UTmp)
  ! TODO: read additional data (e.g. indicators etc), if applicable
ELSE
  ! We need to interpolate the solution to the new computational grid
  SWRITE(UNIT_stdOut,*)'Interpolating base flow from file with N_HDF5=',N_HDF5,' to N=',PP_N
#if PP_dim==3
  N_HDF5Z=N_HDF5
#else
  N_HDF5Z=0
#endif
  ALLOCATE(Vdm_NHDF5_N(0:N_HDF5,0:PP_N))
  ALLOCATE(UTmp2(nVar_HDF5,0:N_HDF5,0:N_HDF5,0:N_HDF5Z,nElems))
  CALL GetVandermonde(N_HDF5,NodeType_HDF5,PP_N,NodeType,Vdm_NHDF5_N,modal=.TRUE.)
  CALL ReadArray('Mean',5,(/nVar_HDF5,N_HDF5+1,N_HDF5+1,N_HDF5Z+1,nElems/),OffsetElem,5,RealArray=UTmp2)
  DO iElem=1,nElems
    CALL ChangeBasis3D(nVar_HDF5,N_HDF5,PP_N,Vdm_NHDF5_N,UTmp2(:,:,:,:,iElem),UTmp(:,:,:,:,iElem))
  END DO
  DEALLOCATE(UTmp2,Vdm_NHDF5_N)
END IF

! Convert to primitive if necessary
IF(primVars)THEN
  UBase(mapPrim(:),:,:,:,:)=UTmp
ELSE 
  KappaM1=Kappa-1.
  UBase(DENS,:,:,:,:)=UTmp(mapCons(1),:,:,:,:)
  UBase(VEL1,:,:,:,:)=UTmp(mapCons(2),:,:,:,:)/UBase(1,:,:,:,:)
  UBase(VEL2,:,:,:,:)=UTmp(mapCons(3),:,:,:,:)/UBase(1,:,:,:,:)
  UBase(VEL3,:,:,:,:)=UTmp(mapCons(4),:,:,:,:)/UBase(1,:,:,:,:)
  UBase(PRES,:,:,:,:)=KappaM1*(UTmp(mapCons(5),:,:,:,:)-0.5*SUM(UBase(VELV,:,:,:,:)*UTmp(mapCons(2:4),:,:,:,:),DIM=1))
END IF
UBase(6,:,:,:,:)=SQRT(Kappa*UBase(5,:,:,:,:)/UBase(1,:,:,:,:))

DEALLOCATE(UTmp,VarNames_HDF5)
CALL CloseDataFile()
SWRITE(UNIT_stdOut,*)'DONE READING BASE FLOW!'
END SUBROUTINE ReadBaseFlow


!==================================================================================================================================
!> Compute gradients and surface contributions of inhomogeneous baseflow
!==================================================================================================================================
SUBROUTINE CalcBaseGradandSurf()
! MODULES
USE MOD_PreProc
USE MOD_Baseflow_Vars,       ONLY:UBase,gradUbx,gradUby,gradUbz,UBase_Master,UBase_Slave
USE MOD_Mesh_Vars,           ONLY:nElems
USE MOD_Mesh_Vars,           ONLY:Metrics_fTilde,Metrics_gTilde,Metrics_hTilde   ! metrics
USE MOD_ProlongToFace,       ONLY:ProlongToFace
USE MOD_Interpolation_Vars,  ONLY: L_Minus,L_Plus,xGP
USE MOD_Basis,               ONLY: PolynomialDerivativeMatrix
USE MOD_FillMortarBF        ,ONLY: U_MortarBF
#if USE_MPI
USE MOD_MPI_Vars
USE MOD_MPI,                 ONLY: StartReceiveMPIData,StartSendMPIData,FinishExchangeMPIData
USE MOD_Mesh_Vars,           ONLY:nSides
#endif /*USE_MPI*/
USE MOD_ApplyJacobianCons, ONLY:ApplyJacobianCons
USE MOD_Mesh_Vars,           ONLY:sJ
USE MOD_Globals ,            ONLY:  myRank
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER   :: iElem,i,j,k,l
REAL,DIMENSION(PP_nVarBase)   :: gradUxi,gradUeta,gradUzeta
REAL      :: D(0:PP_N,0:PP_N),D_T(0:PP_N,0:PP_N)
#if USE_MPI
INTEGER   :: DataSizeSideBF
#endif /*USE_MPI*/
!==================================================================================================================================
! build D matrix
CALL PolynomialDerivativeMatrix(PP_N,xGP,D)
D_T=TRANSPOSE(D)
! Fill gradients (local polynomial derivatives)
DO iElem=1,nElems
  DO k=0,PP_NZ; DO j=0,PP_N; DO i=0,PP_N
    gradUxi     =             D_T(0,i)*UBase(:,0,j,k,iElem)
    gradUeta    =             D_T(0,j)*UBase(:,i,0,k,iElem)
#if (PP_dim==3)
    gradUzeta   =             D_T(0,k)*UBase(:,i,j,0,iElem)
#endif

    DO l=1,PP_N
      gradUxi   = gradUxi   + D_T(l,i)*UBase(:,l,j,k,iElem)
      gradUeta  = gradUeta  + D_T(l,j)*UBase(:,i,l,k,iElem)
#if (PP_dim==3)
      gradUzeta = gradUzeta + D_T(l,k)*UBase(:,i,j,l,iElem)
#endif
    END DO
#if (PP_dim==3)
    gradUBx(:,i,j,k,iElem) = (Metrics_fTilde(1,i,j,k,iElem,0)*gradUxi   &
                           + Metrics_gTilde(1,i,j,k,iElem,0)*gradUeta  &
                           + Metrics_hTilde(1,i,j,k,iElem,0)*gradUzeta)*sJ(i,j,k,iElem,0)
    gradUBy(:,i,j,k,iElem) = (Metrics_fTilde(2,i,j,k,iElem,0)*gradUxi   &
                           + Metrics_gTilde(2,i,j,k,iElem,0)*gradUeta  &
                           + Metrics_hTilde(2,i,j,k,iElem,0)*gradUzeta)*sJ(i,j,k,iElem,0)
    gradUBz(:,i,j,k,iElem) = (Metrics_fTilde(3,i,j,k,iElem,0)*gradUxi   &
                           + Metrics_gTilde(3,i,j,k,iElem,0)*gradUeta  &
                           + Metrics_hTilde(3,i,j,k,iElem,0)*gradUzeta)*sJ(i,j,k,iElem,0)
#else
    gradUBx(:,i,j,k,iElem) = (Metrics_fTilde(1,i,j,k,iElem,0)*gradUxi   &
                           + Metrics_gTilde(1,i,j,k,iElem,0)*gradUeta)*sJ(i,j,k,iElem,0)
    gradUBy(:,i,j,k,iElem) = (Metrics_fTilde(2,i,j,k,iElem,0)*gradUxi   &
                           + Metrics_gTilde(2,i,j,k,iElem,0)*gradUeta)*sJ(i,j,k,iElem,0)
#endif /*PP_dim==3*/
  END DO; END DO; END DO ! i,j,k
END DO ! iElem

! Fill surface data
#if USE_MPI
!DataSizeSideBF=PP_nVarBase*(PP_N+1)**(PP_dim-1)
DataSizeSideBF=PP_nVarBase*(PP_N+1)*(PP_NZ+1)
CALL StartReceiveMPIData(UBase_slave,DataSizeSideBF,1,nSides,MPIRequest_U(:,SEND),SendID=2) ! Receive MINE / U_slave: slave -> master
CALL ProlongToFace(PP_nVarBase,PP_N,UBase,UBase_master,UBase_slave,L_Minus,L_Plus,doMPISides=.TRUE.)
CALL U_MortarBF(UBase_master,UBase_slave,doMPISides=.TRUE.)
CALL StartSendMPIData(UBase_slave,DataSizeSideBF,1,nSides,MPIRequest_U(:,RECV),SendID=2) ! SEND YOUR / U_slave: slave -> master
#endif /*MPI*/
CALL ProlongToFace(PP_nVarBase,PP_N,UBase,UBase_master,UBase_slave,L_Minus,L_Plus,doMPISides=.FALSE.)
#if USE_MPI
CALL FinishExchangeMPIData(2*nNbProcs,MPIRequest_U)        ! UBase_slave: slave -> master 
#endif /*MPI*/
CALL U_MortarBF(UBase_master,UBase_slave,doMPISides=.FALSE.)

END SUBROUTINE CalcBaseGradandSurf

!==================================================================================================================================
!> Finalizes the baseflow
!==================================================================================================================================
SUBROUTINE FinalizeBaseflow()
! MODULES
USE MOD_Baseflow_Vars
IMPLICIT NONE
!==================================================================================================================================
SDEALLOCATE(UBase)
SDEALLOCATE(gradUbx)
SDEALLOCATE(gradUby)
SDEALLOCATE(gradUbz)
SDEALLOCATE(UBase_Master)
SDEALLOCATE(UBase_Slave)
END SUBROUTINE FinalizeBaseflow

END MODULE MOD_Baseflow
