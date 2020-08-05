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
!> This module contains the dependency tables for the quantities to be visualized by the posti routines
!==================================================================================================================================
MODULE MOD_EOS_Posti_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

#if PARABOLIC
#define CUT(x)
INTEGER,PARAMETER :: nVarDepEOS=40
#else
#define CUT(x) x!
INTEGER,PARAMETER :: nVarDepEOS=23
#endif
! ATTENTION: The first     5 variables must be the conservative ones
!            The following 5 variables must be the primitive ones
!           E
!           n
!           e                                    M                                  W
!           r                                    e                                  a
!           g                                    s                                  l
!           y                    E               h                V N               l
!           S            V       n       P       V                o o               F
!           t            e     E t   T   r       e                r r               r W
!           a            l     n h   o   e       l                t m               i a
! W         g            o     e a   t   s       o                i a               c l
! i         n            c V   r l   a T s M M M c                c l         W W W t l
! t         a            i e   g p   l o u e e e i                i i         a a a i H
! h         t         T  t l   y y   T t r s s s t                t z         l l l o e
! D         i         e  y o   S S   e a e h h h y          V V V y e   D Q   l l l n a
! G   M M M o V V V   m  M c   t t   m l T V V V M          o o o M d   i C S F F F M t
! O   o o o n e e e P p  a i   a a   p P i e e e a          r r r a H   l r c r r r a T
! p D m m m D l l l r e  g t   g g E e r m l l l g          t t t g e L a i h i i i g r
! e e e e e e o o o e r  n y   n n n r e e o o o n          i i i n l a t t l c c c n a
! r n n n n n c c c s a  i S   a a t a s D c c c i          c c c i i m a e i t t t i n
! a s t t t s i i i s t  t o M t t r t s e i i i t          i i i t c b t r e i i i t s
! t i u u u i t t t u u  u u a i i o u u r t t t u          t t t u i d i i r o o o u f
! o t m m m t y y y r r  d n c o o p r r i y y y d          y y y d t a o o e n n n d e x y z
! r y X Y Z y X Y Z e e  e d h n n y e e v X Y Z e          X Y Z e y 2 n n n X Y Z e r + + +
INTEGER,DIMENSION(1:nVarDepEOS,0:nVarDepEOS),PARAMETER :: DepTableEOS = TRANSPOSE(RESHAPE(&
(/&
  0,1,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1  Density
  0,0,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2  MomentumX
  0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3  MomentumY
  0,0,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4  MomentumZ
  0,0,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !5  EnergyStagnationDensity
  0,1,1,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !6  VelocityX
  0,1,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !7  VelocityY
  0,1,0,0,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !8  VelocityZ
  0,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !9  Pressure
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !10 Temperature
  0,1,1,1,1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !11 VelocityMagnitude
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !12 VelocitySound
  0,0,0,0,0,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !13 Mach
  0,1,0,0,0,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !14 EnergyStagnation
  0,1,0,0,0,0,0,0,0,1,0, 0,0,0,1,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !15 EnthalpyStagnation
  0,1,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !16 Entropy
  0,0,0,0,0,0,0,0,0,0,1, 1,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !17 TotalTemperature
  0,1,0,0,0,0,0,0,0,1,0, 1,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !18 TotalPressure
  1,1,1,1,1,1,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !19 PressureTimeDeriv
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,1,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !20 MeshVelocityX
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,1,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !21 MeshVelocityY
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,1,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !22 MeshVelocityZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,1  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !23 MeshVelocityMagnitude
#if PARABOLIC
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !24 VorticityX
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !25 VorticityY
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !26 VorticityZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !27 VorticityMagnitude
  1,1,1,1,1,0,0,0,0,0,0, 1,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !28 NormalizedHelicity
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !29 Lambda2
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !30 Dilatation
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !31 QCriterion
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !32 Schlieren
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !33 WallFrictionX
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !34 WallFrictionY
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !35 WallFrictionZ
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0 ,& !36 WallFrictionMagnitude
  1,0,0,0,0,0,0,0,0,0,1, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !37 WallHeatTransfer
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0 ,& !38 x+
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0 ,& !39 y+
  1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0         ,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0  & !40 z+
#endif
/),(/nVarDepEOS+1,nVarDepEOS/)))

! Mark all quantities that can be calculated exclusively on the surface
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepSurfaceOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 &
/)

! Mark all quantities that can be calculated exclusively in the volume and must be prolonged to the surface from the volume
INTEGER,DIMENSION(1:nVarDepEOS),PARAMETER :: DepVolumeOnlyEOS = &
(/  0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,1,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 &
/)

#if FV_ENABLED && FV_RECONSTRUCT
!           E
!           n
!           e
!           r
!           g
!           y
!           S
!           t
!           a
! W         g
! i         n
! t         a
! h         t         T
! D         i         e
! G   M M M o V V V   m
! O   o o o n e e e P p
! p D m m m D l l l r e
! e e e e e e o o o e r
! r n n n n n c c c s a
! a s t t t s i i i s t
! t i u u u i t t t u u
! o t m m m t y y y r r
! r y X Y Z y X Y Z e e
INTEGER,DIMENSION(PP_nVar,0:nVarDepEOS),PARAMETER :: DepTablePrimToCons =TRANSPOSE(RESHAPE(&
(/&
  0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !1 Density
  0,1,0,0,0,0,1,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !2 MomentumX
  0,1,0,0,0,0,0,1,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !3 MomentumY
  0,1,0,0,0,0,0,0,1,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0, CUT(&)  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ,& !4 MomentumZ
  0,1,0,0,0,0,1,1,1,1,0, 0,0,0,0,0,0,0,0,0,0,0,0,0  CUT(&) ,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0  & !5 EnergyStagnationDensity
/),(/nVarDepEOS+1,5/)))
#endif
#undef CUT

CHARACTER(LEN=255),DIMENSION(nVarDepEOS),PARAMETER :: DepNames = &
(/ CHARACTER(LEN=255) ::    &
"Density"                  ,& !1
"MomentumX"                ,& !2
"MomentumY"                ,& !3
"MomentumZ"                ,& !4
"EnergyStagnationDensity"  ,& !5
"VelocityX"                ,& !6
"VelocityY"                ,& !7
"VelocityZ"                ,& !8
"Pressure"                 ,& !9
"Temperature"              ,& !10
"VelocityMagnitude"        ,& !11
"VelocitySound"            ,& !12
"Mach"                     ,& !13
"EnergyStagnation"         ,& !14
"EnthalpyStagnation"       ,& !15
"Entropy"                  ,& !16
"TotalTemperature"         ,& !17
"TotalPressure"            ,& !18
"PressureTimeDeriv"        ,& !19
"MeshVelocityX"            ,& !20
"MeshVelocityY"            ,& !21
"MeshVelocityZ"            ,& !22
"MeshVelocityMagnitude"     & !23
#if PARABOLIC
,"VorticityX"              ,& !24
"VorticityY"               ,& !25
"VorticityZ"               ,& !26
"VorticityMagnitude"       ,& !27
"NormalizedHelicity"       ,& !28
"Lambda2"                  ,& !29
"Dilatation"               ,& !30
"QCriterion"               ,& !31
"Schlieren"                ,& !32
"WallFrictionX"            ,& !33
"WallFrictionY"            ,& !34
"WallFrictionZ"            ,& !35
"WallFrictionMagnitude"    ,& !36
"WallHeatTransfer"         ,& !37
"x+"                       ,& !38
"y+"                       ,& !39
"z+"                        & !40
#endif /*PARABOLIC*/
/)

END MODULE MOD_EOS_Posti_Vars
