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

!==================================================================================================================================
!> Contains global variables provided by the movemesh routines.
!==================================================================================================================================
MODULE MOD_MoveMesh_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
! User defined parameters
INTEGER                        :: MeshMoveType                   !< Type of mesh movement
! Mappings
INTEGER,ALLOCATABLE            :: S2V2_NGeo(:,:,:,:,:)           !< Side to volume mapping on NGeo
INTEGER,ALLOCATABLE            :: FS2M_NGeo(:,:,:,:)             !< Flip side to mortar mapping on NGeo
! Arrays containing mesh velocites
REAL,ALLOCATABLE               :: v_NGeo(:,:,:,:,:)              !< uvw velocities (first index 1:3) of the volume mesh points
REAL,ALLOCATABLE               :: v_NGeo_Face(:,:,:,:)           !< uvw velocities (first index 1:3) of the surface mesh points
REAL,ALLOCATABLE               :: Elem_vGP(:,:,:,:,:)            !< uvw velocities (first index 1:3) of the volume Gauss Point
REAL,ALLOCATABLE               :: Face_vGP(:,:,:,:)              !< uvw velocities (first index 1:3) of the surface Gauss Point
! Arrays needed for time dependent mesh representation
REAL,ALLOCATABLE               :: NodeCoords_actual(:,:,:,:,:)   !< Current Node Coordinates on CL points (NGeo)
REAL,ALLOCATABLE               :: Vdm_CLNgeo_GaussN(:,:)         !< Vandermonde from CL points on Ngeo to Gauss points on N
! Variables specific to mesh move types
LOGICAL                        :: MeshVelAlreadySet = .FALSE.    !< Switch to set constant velocities only once
LOGICAL,ALLOCATABLE            :: DoMeshAccelerate(:)            !< Flag to allow for movement with constant acceleration
REAL,ALLOCATABLE               :: MeshAcc(:,:)                   !< Mesh velocity array
REAL,ALLOCATABLE               :: MeshVel(:,:)                   !< Mesh velocity array
REAL,ALLOCATABLE               :: RotationAngVel(:)              !< Angular velocity of rotation
REAL,ALLOCATABLE               :: RotationCenter(:,:)            !< Center of rotational mesh movement
#if FV_ENABLED
! FV specific variables
REAL,ALLOCATABLE               :: Vdm_CLNGeo_FVBdryx(:,:)
REAL,ALLOCATABLE               :: Vdm_CLNGeo_FV     (:,:)
#endif
! Mortar specific
INTEGER                        :: nSmallInnerMortarSides         !< Number of inner small mortar sides
INTEGER,ALLOCATABLE,TARGET     :: SmallInnerMortarSideIDs(:)     !< List of inner small mortar side IDs
! Send operation (from my big mortars)
INTEGER                        :: nSmallMortarMPISidesSend       !< Number of small mortar sides where we send informations
INTEGER                        :: nNbProcsMortarSend             !< Number of neighbour procs
INTEGER,ALLOCATABLE            :: NbProcsMortarSend(:)           !< Ranks of neighbour procs
INTEGER,ALLOCATABLE            :: nMortarSidesSend(:)            !< Number of sides to send per rank
INTEGER,ALLOCATABLE            :: offsetMortarSidesSend(:)       !< Offset in send array
INTEGER,ALLOCATABLE            :: mapMortarSidesSend(:)          !< Map from entry in send array to real side IDs
! Recieve operation (to my small mortars)
INTEGER                        :: nSmallMortarMPISidesRcv        !< Number of small mortar where we need to recieve informations
INTEGER,ALLOCATABLE,TARGET     :: SmallMortarMPISidesRcv(:,:)    !< List of small send side IDs
INTEGER                        :: nNbProcsMortarRcv              !< Number of neighbour procs
INTEGER,ALLOCATABLE            :: NbProcsMortarRcv(:)            !< Ranks of neighbour procs
INTEGER,ALLOCATABLE            :: nMortarSidesRcv(:)             !< Number of sides to recive per rank
INTEGER,ALLOCATABLE            :: offsetMortarSidesRcv(:)        !< Offset in recieve array
INTEGER,ALLOCATABLE            :: mapMortarSidesRcv(:)           !< Map from entry in recieve array to real side IDs

LOGICAL                        :: MeshMoveInitIsDone =.FALSE.    !< Switch to check execution of init function
END MODULE MOD_MoveMesh_Vars
