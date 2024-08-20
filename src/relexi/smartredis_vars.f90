!=================================================================================================================================
! Copyright (c) 2010-2022 Prof. Claus-Dieter Munz
! Copyright (c) 2022-2024 Prof. Andrea Beck
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
!==================================================================================================================================
!> Contains the variables of necessary for the SmartRedis adapter that should be globally accessible!
!==================================================================================================================================
MODULE MOD_SmartRedis_Vars
#if USE_SMARTREDIS
! MODULES
USE SmartRedis_Client,  ONLY: client_type
IMPLICIT NONE
PUBLIC
SAVE

!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
TYPE(client_type)   :: client           ! Client instance of SmartRedis to communicate with Redis Database
LOGICAL             :: dbIsClustered    ! Indicate whether the Redis Database is clustered, i.e. distributed on different nodes
LOGICAL             :: doSmartRedis     ! Flag whether communication with SmartRedis should be done
LOGICAL             :: useInvariants    ! If true, the invariants of the gradient tensor are used as state instead of the velocities
LOGICAL             :: doNormInvariants ! Normalize invariants by first one
INTEGER             :: SR_nVarAction    ! Number/Dimension of actions received by agent per element
INTEGER             :: SR_Error         ! Integer containing the SmartRedis Error
#endif

END MODULE MOD_SmartRedis_Vars
