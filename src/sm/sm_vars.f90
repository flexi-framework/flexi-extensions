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

!==================================================================================================================================
!> Variables used for sliding mesh
!==================================================================================================================================
MODULE MOD_SM_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                    :: SMInitIsDone=.FALSE.             !> 
REAL,ALLOCATABLE,TARGET    :: M_SM_0_1(:,:,:),M_SM_0_2(:,:,:)  !> 
REAL,ALLOCATABLE,TARGET    :: M_SM_1_0(:,:,:),M_SM_2_0(:,:,:)  !>   
INTEGER,ALLOCATABLE        :: SMFlip(:)                        !>
INTEGER,ALLOCATABLE        :: SM_nMPIMortars_Proc(:)           !> 
INTEGER,ALLOCATABLE,TARGET :: SMMortarToRotSide(:,:)           !> 
INTEGER,ALLOCATABLE,TARGET :: SMMortarToStatSide(:,:)          !> 
REAL,ALLOCATABLE           :: NormVec_SM (:,:,:,:)             !> 
REAL,ALLOCATABLE           :: TangVec1_SM(:,:,:,:)             !> 
REAL,ALLOCATABLE           :: TangVec2_SM(:,:,:,:)             !> 
REAL,ALLOCATABLE           :: U_MorStat(:,:,:,:)               !> 
REAL,ALLOCATABLE           :: U_MorRot(:,:,:,:)                !>
REAL,ALLOCATABLE           :: Flux_MorStat(:,:,:,:)            !> 
REAL,ALLOCATABLE           :: Flux_MorRot(:,:,:,:)             !> 
REAL,ALLOCATABLE           :: UPrim_MorStat(:,:,:,:)           !>
REAL,ALLOCATABLE           :: UPrim_MorRot (:,:,:,:)           !>
#if PARABOLIC
REAL,ALLOCATABLE           :: gradUx_MorStat(:,:,:,:)          !>
REAL,ALLOCATABLE           :: gradUx_MorRot (:,:,:,:)          !>
REAL,ALLOCATABLE           :: gradUy_MorStat(:,:,:,:)          !>
REAL,ALLOCATABLE           :: gradUy_MorRot (:,:,:,:)          !>
REAL,ALLOCATABLE           :: gradUz_MorStat(:,:,:,:)          !>
REAL,ALLOCATABLE           :: gradUz_MorRot (:,:,:,:)          !>
#if EDDYVISCOSITY
REAL,ALLOCATABLE           :: muSGS_MorStat(:,:,:,:)           !>
REAL,ALLOCATABLE           :: muSGS_MorRot( :,:,:,:)           !>
#endif /*EDDYVISCOSITY*/
#endif /*PARABOLIC*/
#if FV_ENABLED
INTEGER,ALLOCATABLE        :: SM_Elems(:)                      !>
#endif /* FV_ENABLED */
REAL,ALLOCATABLE           :: Face_vGPSM(:,:,:)                !>
INTEGER,ALLOCATABLE        :: SMn_prev(:)                      !> 

INTEGER                    :: nSlidingMeshInterfaces
INTEGER                    :: nSlidingMeshPartitions
INTEGER,ALLOCATABLE        :: SlidingMeshType(:)               !>
INTEGER,ALLOCATABLE        :: SlidingMeshDirection(:)          !>
REAL,ALLOCATABLE           :: SlidingMeshBoundaries(:,:)       !>
!==================================================================================================================================
END MODULE MOD_SM_Vars
