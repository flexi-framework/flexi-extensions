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
!> Special module for mortar operation on mesh coordinates. Is in a seperate module since we need different mortar operators.
!===================================================================================================================================
MODULE MOD_FillMortarMesh
IMPLICIT NONE
PRIVATE

INTERFACE Mesh_Mortar
  MODULE PROCEDURE Mesh_Mortar
END INTERFACE

#if GCL
INTERFACE Flux_MortarGCL
  MODULE PROCEDURE Flux_MortarGCL
END INTERFACE
#endif

PUBLIC::Mesh_Mortar
#if GCL
PUBLIC::Flux_MortarGCL
#endif

CONTAINS
!==================================================================================================================================
!> Fills small non-conforming sides with data from the corresponding large side, using 1D interpolation operators M_0_1,M_0_2.
!> This routine is used to interpolate the mesh displacements and velocities.
!> These are continuous over cell interfaces, so we need no master and slave arrays. Also this is always done on the CL points,
!> even for DG/FV or FV/FV interfaces. So no special FV algorithm is needed.
!>
!>       Type 1               Type 2              Type3
!>        eta                  eta                 eta
!>         ^                    ^                   ^
!>         |                    |                   |
!>     +---+---+            +---+---+           +---+---+
!>     | 3 | 4 |            |   2   |           |   |   |
!>     +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>     | 1 | 2 |            |   1   |           |   |   |
!>     +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE Mesh_Mortar(MeshData_in,doMPISides)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars,  ONLY: M_0_1_Mesh,M_0_2_Mesh
USE MOD_Mesh_Vars,    ONLY: MortarType,MortarInfo,NGeo
USE MOD_Mesh_Vars,    ONLY: firstMortarInnerSide,lastMortarInnerSide
USE MOD_Mesh_Vars,    ONLY: firstMortarMPISide,lastMortarMPISide
USE MOD_Mesh_Vars,    ONLY: nSides
USE MOD_MoveMesh_Vars,ONLY: FS2M_NGeo
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: MeshData_in(1:3,0:NGeo,0:ZDIM(NGeo),1:nSides) !< (INOUT) can be node coordinates or velocites
LOGICAL,INTENT(IN) :: doMPISides                                  !< flag whether MPI sides are processed
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,locSide,flip
REAL         :: MeshData_tmp( 1:3,0:NGeo,0:ZDIM(NGeo),1:4)
#if PP_dim == 3
REAL         :: MeshData_tmp2(1:3,0:NGeo,0:ZDIM(NGeo),1:2)
#endif
REAL,POINTER :: M1(:,:),M2(:,:)
!==================================================================================================================================
!                         doMPISides==True   doMPISides==False
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)

M1=>M_0_1_Mesh; M2=>M_0_2_Mesh

DO MortarSideID=firstMortarSideID,lastMortarSideID

#if PP_dim == 3
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    !first in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    MeshData_tmp2(iVar,p,:,1)  =  M1 * MeshData_in(iVar,p,:,MortarSideID)
    !    MeshData_tmp2(iVar,p,:,2)  =  M2 * MeshData_in(iVar,p,:,MortarSideID)
    DO q=0,NGeo
      DO p=0,NGeo ! for every xi-layer perform Mortar operation in eta-direction 
        MeshData_tmp2(:,p,q,1)=                  M1(0,q)*MeshData_in(:,p,0,MortarSideID)
        MeshData_tmp2(:,p,q,2)=                  M2(0,q)*MeshData_in(:,p,0,MortarSideID)
        DO l=1,NGeo
          MeshData_tmp2(:,p,q,1)=MeshData_tmp2(:,p,q,1)+M1(l,q)*MeshData_in(:,p,l,MortarSideID)
          MeshData_tmp2(:,p,q,2)=MeshData_tmp2(:,p,q,2)+M2(l,q)*MeshData_in(:,p,l,MortarSideID)
        END DO
      END DO
    END DO
    ! then in xi
    DO q=0,NGeo ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    MeshData_tmp(iVar,:,q,1)  =  M1 * MeshData_tmp2(iVar,:,q,1)
      !    MeshData_tmp(iVar,:,q,2)  =  M2 * MeshData_tmp2(iVar,:,q,1)
      !    MeshData_tmp(iVar,:,q,3)  =  M1 * MeshData_tmp2(iVar,:,q,2)
      !    MeshData_tmp(iVar,:,q,4)  =  M2 * MeshData_tmp2(iVar,:,q,2)
      DO p=0,NGeo
        MeshData_tmp(:,p,q,1)=                 M1(0,p)*MeshData_tmp2(:,0,q,1)
        MeshData_tmp(:,p,q,2)=                 M2(0,p)*MeshData_tmp2(:,0,q,1)
        MeshData_tmp(:,p,q,3)=                 M1(0,p)*MeshData_tmp2(:,0,q,2)
        MeshData_tmp(:,p,q,4)=                 M2(0,p)*MeshData_tmp2(:,0,q,2)
        DO l=1,NGeo
          MeshData_tmp(:,p,q,1)=MeshData_tmp(:,p,q,1)+M1(l,p)*MeshData_tmp2(:,l,q,1)
          MeshData_tmp(:,p,q,2)=MeshData_tmp(:,p,q,2)+M2(l,p)*MeshData_tmp2(:,l,q,1)
          MeshData_tmp(:,p,q,3)=MeshData_tmp(:,p,q,3)+M1(l,p)*MeshData_tmp2(:,l,q,2)
          MeshData_tmp(:,p,q,4)=MeshData_tmp(:,p,q,4)+M2(l,p)*MeshData_tmp2(:,l,q,2)
        END DO !l=1,NGeo
      END DO
    END DO 

  CASE(2) !1->2 in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    MeshData_tmp(iVar,p,:,1)  =  M1 * MeshData_in(iVar,p,:,MortarSideID)
    !    MeshData_tmp(iVar,p,:,2)  =  M2 * MeshData_in(iVar,p,:,MortarSideID)
    DO q=0,NGeo
      DO p=0,NGeo ! for every xi-layer perform Mortar operation in eta-direction
        MeshData_tmp(:,p,q,1)=                 M1(0,q)*MeshData_in(:,p,0,MortarSideID)
        MeshData_tmp(:,p,q,2)=                 M2(0,q)*MeshData_in(:,p,0,MortarSideID)
        DO l=1,NGeo
          MeshData_tmp(:,p,q,1)=MeshData_tmp(:,p,q,1)+M1(l,q)*MeshData_in(:,p,l,MortarSideID)
          MeshData_tmp(:,p,q,2)=MeshData_tmp(:,p,q,2)+M2(l,q)*MeshData_in(:,p,l,MortarSideID)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi      NOTE: In 2D only the first space index can be a Mortar (second index is always 0)!!!
#endif
    DO q=0,ZDIM(NGeo) ! for every eta-layer perform Mortar operation in xi-direction
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    MeshData_tmp(iVar,:,q,1)  =  M1 * MeshData_in(iVar,:,q,MortarSideID)
      !    MeshData_tmp(iVar,:,q,2)  =  M2 * MeshData_in(iVar,:,q,MortarSideID)
      DO p=0,NGeo
        MeshData_tmp(:,p,q,1)=                 M1(0,p)*MeshData_in(:,0,q,MortarSideID)
        MeshData_tmp(:,p,q,2)=                 M2(0,p)*MeshData_in(:,0,q,MortarSideID)
        DO l=1,NGeo
          MeshData_tmp(:,p,q,1)=MeshData_tmp(:,p,q,1)+M1(l,p)*MeshData_in(:,l,q,MortarSideID)
          MeshData_tmp(:,p,q,2)=MeshData_tmp(:,p,q,2)+M2(l,p)*MeshData_in(:,l,q,MortarSideID)
        END DO
      END DO
    END DO
#if PP_dim == 3
  END SELECT ! mortarType(SideID)
#endif

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  locSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID= MortarInfo(MI_SIDEID,iMortar,locSide)
    flip  = MortarInfo(MI_FLIP,iMortar,locSide)
    SELECT CASE(flip)
      CASE(0) ! master side
        MeshData_in(:,:,:,SideID)=MeshData_tmp(:,:,:,iMortar)
      CASE(1:4) ! slave side
        DO q=0,ZDIM(NGeo); DO p=0,NGeo
          MeshData_in(:,p,q,SideID)=MeshData_tmp(:,FS2M_NGeo(1,p,q,flip), &
                                           FS2M_NGeo(2,p,q,flip),iMortar)
        END DO; END DO ! q, p
    END SELECT !flip(iMortar)
  END DO !iMortar
END DO !MortarSideID
END SUBROUTINE Mesh_Mortar

#if GCL
!==================================================================================================================================
!>  Fills master side from small non-conforming sides, using 1D projection operators M_1_0,M_2_0
!>
!> This routine is used to project the numerical flux at the small sides of the nonconforming interface to the corresponding large
!>  ones.
!>
!>        Type 1               Type 2              Type3
!>         eta                  eta                 eta
!>          ^                    ^                   ^
!>          |                    |                   |
!>      +---+---+            +---+---+           +---+---+
!>      | 3 | 4 |            |   2   |           |   |   |
!>      +---+---+ --->  xi   +---+---+ --->  xi  + 1 + 2 + --->  xi
!>      | 1 | 2 |            |   1   |           |   |   |
!>      +---+---+            +---+---+           +---+---+
!>
!==================================================================================================================================
SUBROUTINE Flux_MortarGCL(Flux_master,Flux_slave,doMPISides,weak,onlyFV)
! MODULES
USE MOD_Preproc
USE MOD_Mortar_Vars, ONLY: M_1_0,M_2_0
#if FV_ENABLED
USE MOD_Mortar_Vars, ONLY: FV_M_1_0,FV_M_2_0
#endif
USE MOD_Mesh_Vars,   ONLY: MortarType,MortarInfo,nSides
USE MOD_Mesh_Vars,   ONLY: firstMortarInnerSide,lastMortarInnerSide,FS2M
USE MOD_Mesh_Vars,   ONLY: firstMortarMPISide,lastMortarMPISide
#if FV_ENABLED
USE MOD_FV_Vars     ,ONLY: FV_Elems_master
#endif
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL,INTENT(INOUT) :: Flux_master(1:1,0:PP_N,0:PP_NZ,1:nSides) !< input data
REAL,INTENT(IN)    :: Flux_slave (1:1,0:PP_N,0:PP_NZ,1:nSides) !< input data
LOGICAL,INTENT(IN) :: doMPISides                             !< flag whether MPI sides are processed
LOGICAL,INTENT(IN) :: weak                                   !< flag whether strong or weak form is used
LOGICAL,INTENT(IN),OPTIONAL :: onlyFV                        !< flag whether only FV sides should be handled
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER      :: p,q,l
INTEGER      :: iMortar,nMortars
INTEGER      :: firstMortarSideID,lastMortarSideID
INTEGER      :: MortarSideID,SideID,iSide,flip
REAL         :: Flux_tmp( 1,0:PP_N,0:PP_NZ,1:4)
#if PP_dim == 3
REAL         :: Flux_tmp2(1,0:PP_N,0:PP_NZ,1:2)
#endif
REAL,POINTER :: M1(:,:),M2(:,:)
LOGICAL      :: onlyFV_loc
!==================================================================================================================================
!                         doMPISides==True   doMPISides==False
firstMortarSideID = MERGE(firstMortarMPISide,firstMortarInnerSide,doMPISides)
 lastMortarSideID = MERGE( lastMortarMPISide, lastMortarInnerSide,doMPISides)
onlyFV_loc = MERGE(onlyFV, .FALSE., PRESENT(onlyFV))

M1=>M_1_0; M2=>M_2_0
DO MortarSideID=firstMortarSideID,lastMortarSideID
#if FV_ENABLED
  IF (onlyFV_loc) THEN
    IF(FV_Elems_master(MortarSideID).EQ.0) CYCLE ! DG
  END IF
#endif

  nMortars=MERGE(4,2,MortarType(1,MortarSideID).EQ.1)
  iSide=MortarType(2,MortarSideID)
  DO iMortar=1,nMortars
    SideID = MortarInfo(MI_SIDEID,iMortar,iSide)
    flip   = MortarInfo(MI_FLIP,iMortar,iSide)
    SELECT CASE(flip)
    CASE(0) ! master side
      Flux_tmp(:,:,:,iMortar)=Flux_master(:,:,:,SideID)
    CASE(1:4) ! slave sides (should only occur for MPI)
      IF(weak)THEN
        DO q=0,PP_NZ; DO p=0,PP_N
          Flux_tmp(:,p,q,iMortar)=-Flux_slave(:,FS2M(1,p,q,flip),FS2M(2,p,q,flip),SideID)
        END DO; END DO
      ELSE    ! do not change sign if strong form
        DO q=0,PP_NZ; DO p=0,PP_N
          Flux_tmp(:,p,q,iMortar)= Flux_slave(:,FS2M(1,p,q,flip),FS2M(2,p,q,flip),SideID)
        END DO; END DO
      END IF
    END SELECT
  END DO

#if FV_ENABLED
  IF(FV_Elems_master(MortarSideID).GT.0)THEN ! FV element
    M1=>FV_M_1_0
    M2=>FV_M_2_0
  ELSE ! DG element
    M1=>   M_1_0
    M2=>   M_2_0
  END IF
#endif
#if PP_dim == 3
  SELECT CASE(MortarType(1,MortarSideID))
  CASE(1) !1->4
    ! first in xi
    DO q=0,PP_N ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are four MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_tmp2(iVar,:,q,1)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      !    Flux_tmp2(iVar,:,q,2)  =  M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,PP_N
        Flux_tmp2(:,p,q,1)=                       M1(0,p)*Flux_tmp(:,0,q,1)+M2(0,p)*Flux_tmp(:,0,q,2)
        Flux_tmp2(:,p,q,2)=                       M1(0,p)*Flux_tmp(:,0,q,3)+M2(0,p)*Flux_tmp(:,0,q,4)
        DO l=1,PP_N
          Flux_tmp2(:,p,q,1)=Flux_tmp2(:,p,q,1) + M1(l,p)*Flux_tmp(:,l,q,1)+M2(l,p)*Flux_tmp(:,l,q,2)
          Flux_tmp2(:,p,q,2)=Flux_tmp2(:,p,q,2) + M1(l,p)*Flux_tmp(:,l,q,3)+M2(l,p)*Flux_tmp(:,l,q,4)
        END DO
      END DO
    END DO
    !then in eta
    ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
    !    Flux_master(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp2(iVar,p,:,1) + M2 * Flux_tmp2(iVar,p,:,2)
    DO q=0,PP_N
      DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        Flux_master(:,p,q,MortarSideID)=                                    M1(0,q)*Flux_tmp2(:,p,0,1)+M2(0,q)*Flux_tmp2(:,p,0,2)
        DO l=1,PP_N
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) + M1(l,q)*Flux_tmp2(:,p,l,1)+M2(l,q)*Flux_tmp2(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(2) !1->2 in eta
    ! TODO why not q-loop first?
    DO p=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
      ! The following q- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_master(iVar,p,:,MortarSideID)  =  M1 * Flux_tmp(iVar,p,:,1) + M2 * Flux_tmp(iVar,p,:,2)
      DO q=0,PP_N ! for every xi-layer perform Mortar operation in eta-direction 
        Flux_master(:,p,q,MortarSideID)=                                    M1(0,q)*Flux_tmp(:,p,0,1)+M2(0,q)*Flux_tmp(:,p,0,2)
        DO l=1,PP_N
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) + M1(l,q)*Flux_tmp(:,p,l,1)+M2(l,q)*Flux_tmp(:,p,l,2)
        END DO
      END DO
    END DO

  CASE(3) !1->2 in xi   NOTE: In 2D only XI mortars
#endif
    DO q=0,PP_NZ ! for every eta-layer perform Mortar operation in xi-direction 
      ! The following p- and l-loop are two MATMULs: (ATTENTION M1 and M2 are already transposed in mortar.f90)
      !    Flux_master(iVar,:,q,MortarSideID)  =   M1 * Flux_tmp(iVar,:,q,1) + M2 * Flux_tmp(iVar,:,q,2)
      DO p=0,PP_N
        Flux_master(:,p,q,MortarSideID)=                                    M1(0,p)*Flux_tmp(:,0,q,1)+M2(0,p)*Flux_tmp(:,0,q,2)
        DO l=1,PP_N
          Flux_master(:,p,q,MortarSideID)=Flux_master(:,p,q,MortarSideID) + M1(l,p)*Flux_tmp(:,l,q,1)+M2(l,p)*Flux_tmp(:,l,q,2)
        END DO
      END DO
    END DO

#if PP_dim == 3
  END SELECT ! mortarType(MortarSideID)
#endif
END DO !MortarSideID
END SUBROUTINE Flux_MortarGCL
#endif /*GCL*/

END MODULE MOD_FillMortarMesh
