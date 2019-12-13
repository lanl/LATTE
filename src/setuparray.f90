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

MODULE SETUPARRAY

  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER, ALLOCATABLE :: ELEMPOINTER(:)
  INTEGER, ALLOCATABLE :: MATINDLIST(:), SPININDLIST(:)
  INTEGER, ALLOCATABLE :: BTYPE_INT(:,:)
  INTEGER, ALLOCATABLE :: ORBITAL_LIST(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: CR(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HR0(:)
  REAL(LATTEPREC), ALLOCATABLE :: HES(:), HEP(:), HED(:), HEF(:), ATOCC(:)
  REAL(LATTEPREC), ALLOCATABLE :: H(:,:), BO(:,:), BOZERO(:), H0(:,:), HDIAG(:)
  REAL(LATTEPREC), ALLOCATABLE :: H_ONSITE(:)
  REAL(LATTEPREC), ALLOCATABLE :: ORTHORHO(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: F(:,:), FPP(:,:), FTOT(:,:), FCOUL(:,:), FPLUSD(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: FPUL(:,:), FSCOUL(:,:), FSSPIN(:,:), FSLCN(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: DELTAQ(:), MYCHARGE(:), QLIST(:), OLDQLIST(:)
  REAL(LATTEPREC), ALLOCATABLE :: LCNSHIFT(:)
  REAL(LATTEPREC), ALLOCATABLE :: HUBBARDU(:)
  REAL(LATTEPREC), ALLOCATABLE :: RESPCHI(:)
  REAL(LATTEPREC), ALLOCATABLE :: CUTOFF_LIST(:,:)
  CHARACTER(LEN=1), ALLOCATABLE :: RELAXATOM(:,:)
  CHARACTER(LEN=2), ALLOCATABLE :: ELE(:), ELE1(:), ELE2(:), ATELE(:)
  CHARACTER(LEN=3), ALLOCATABLE :: BTYPE(:)
  CHARACTER(LEN=4), ALLOCATABLE :: BASIS(:)

  ! More Coulomb related data

  REAL(LATTEPREC), ALLOCATABLE :: COULOMBV(:)

END MODULE SETUPARRAY
