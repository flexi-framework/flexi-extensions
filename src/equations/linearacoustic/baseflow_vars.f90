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
!> Contains all variables associated with the mean or base flow of the linear acoustic equation
!==================================================================================================================================
MODULE MOD_Baseflow_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                  :: varMeanFlow              !< variable mean flow: source term needs to be included
INTEGER                  :: BaseflowType
REAL                     :: BaseState(PP_nVarBase)   !< base state in primitive variables (rho,u,v,w,p)
REAL                     :: uShear                   !< shear-layer peak velocity analytical baseflow 2
REAL                     :: dShear                   !< shear-layer thickness analytical baseflow 2
REAL                     :: u1
REAL                     :: u2
REAL,ALLOCATABLE,TARGET  :: UBase(:,:,:,:,:)        
REAL,ALLOCATABLE         :: gradUbx(:,:,:,:,:)
REAL,ALLOCATABLE         :: gradUby(:,:,:,:,:)
REAL,ALLOCATABLE         :: gradUbz(:,:,:,:,:)
REAL,ALLOCATABLE         :: UBase_Master(:,:,:,:)
REAL,ALLOCATABLE         :: UBase_Slave(:,:,:,:)
!==================================================================================================================================
END MODULE MOD_Baseflow_Vars
