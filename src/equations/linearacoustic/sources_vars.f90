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
MODULE MOD_AcSources_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------
INTEGER                        :: AcSourceVariable         !< pressuretimederiv (0), lambvector (1), pressurematderiv (2)
CHARACTER(LEN=255)             :: AcSourceProjectName      !< Project name of source database
CHARACTER(LEN=255)             :: AcSourcePath             !< Source path indicating the source data base to be used
CHARACTER(LEN=255),ALLOCATABLE :: SrcFileNames(:)          !< Source file list
REAL,ALLOCATABLE               :: Time_acsrc(:)            !< Time of complete source file list

CHARACTER(LEN=255),DIMENSION(3),PARAMETER :: VarNamesSrcLamb =&
  (/ CHARACTER(LEN=255) :: 'LambVecX','LambVecY','LambVecZ'/) !< source variable names Lamb vector
CHARACTER(LEN=255),DIMENSION(1),PARAMETER :: VarNamesSrcPTime =&
  (/ CHARACTER(LEN=255) :: 'PressureTimeDeriv'/) !< source variable names pressure time derivative
CHARACTER(LEN=255),DIMENSION(1),PARAMETER :: VarNamesSrcPMat  =&
  (/ CHARACTER(LEN=255) :: 'PressureMatDeriv'/) !< source variable names pressure material derivative

CHARACTER(LEN=255),ALLOCATABLE :: VarNamesSrc(:)
INTEGER                        :: nVarSrc,nVarSrcRead,nVarSrc_HDF5,nSrcFiles
INTEGER                        :: NTime
INTEGER                        :: nextGlob,nextBuffer
INTEGER,ALLOCATABLE            :: varMapRead(:),inputMapRead(:),mapSrcToVar(:)
INTEGER,ALLOCATABLE            :: iBuffer(:),iGlob(:)
REAL,ALLOCATABLE               :: acSrcBuffer(:,:,:,:,:,:),acSrcMean(:,:,:,:,:)
REAL,ALLOCATABLE               :: wBaryAcSrc(:)
INTEGER,DIMENSION(4,4)         :: Stencils
LOGICAL                        :: doMaskSource,doTemporalRamp
REAL                           :: widthMask,rMask
REAL                           :: t0p5,t1
REAL                           :: TEndAc
LOGICAL          :: AcMaskViz             !< Turn on to write a visualization file of the AcMAsk region and strength
INTEGER          :: nAcMaskElems          !< number of elements for which AcMask is applied
INTEGER,ALLOCATABLE :: AcMaskMap(:)       !< mapping from Elem -> AcMaskElem
REAL,ALLOCATABLE :: AcMaskMat(:,:,:,:)    !< precomputed sponge functions per DOF and AcMAsk elem
REAL,ALLOCATABLE,TARGET :: AcMaskBaseFlow(:,:,:,:,:) !< precompute global reference state for whole field

LOGICAL                        :: SkipElemByPoly
INTEGER                        :: nSkipZones
INTEGER,ALLOCATABLE            :: nVertices(:)
REAL,ALLOCATABLE               :: SkipZoneVertex(:,:,:)
INTEGER                        :: nSkippedElems
INTEGER,ALLOCATABLE            :: skippedElems(:)
!==================================================================================================================================
END MODULE MOD_AcSources_Vars
