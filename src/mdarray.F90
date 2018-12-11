!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MDARRAY

  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER :: ISET
  REAL(LATTEPREC), ALLOCATABLE :: V(:,:), MASS(:)
  REAL(LATTEPREC), ALLOCATABLE :: THIST(:), PHIST(:), EHIST(:), VHIST(:)
  REAL(LATTEPREC), ALLOCATABLE :: PHISTX(:), PHISTY(:), PHISTZ(:)
  REAL(LATTEPREC) :: DT, TTARGET, PTARGET, TEMPERATURE, DTZERO, TZERO
  REAL(LATTEPREC) :: AVET, HG, P0, V0, E0
  INTEGER :: CONTITER, NVTON, NPTON, AVEPER, TOINITTEMP, GETHUG
  INTEGER :: MAXITER, DUMPFREQ, RSFREQ, WRTFREQ, THERMPER, THERMRUN
  INTEGER :: ENTROPYITER
  CHARACTER(LEN=10) :: RNDIST, SEEDINIT, NPTTYPE

  ! These are for the Hugoniostat

  INTEGER :: SHOCKON, SHOCKSTART, SHOCKSTOP, SHOCKDIR 
  REAL(LATTEPREC) :: UPARTICLE, USHOCK, C0

!!$  Thermostating
  INTEGER :: SEEDTH
  INTEGER :: SETTH
  REAL(LATTEPREC) :: FRICTION

!!$  These are for Langevin dynamics

  REAL(LATTEPREC), ALLOCATABLE :: FRANPREV(:,:)

!!$  These are for Andersen thermostat

  REAL(LATTEPREC) :: CUMDT
  !  REAL(LATTEPREC) :: CUMDT, TAU
  !  INTEGER :: QITERAND, QITERIN

!!$  These are for Nose thermostat

  REAL(LATTEPREC) :: GAMMA, DGAMMA
  REAL(LATTEPREC) :: CONSMOT

END MODULE MDARRAY
