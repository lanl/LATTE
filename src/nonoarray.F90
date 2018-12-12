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

MODULE NONOARRAY

  USE MYPRECISION 

  IMPLICIT NONE
  SAVE

  ! Work arrays are for DSYEV when generating s^-1/2
  ! SMAT contains the overlap matrix S
  ! UMAT contains eigenvectors of S and NONO_EVALS its eigenvalues
  ! XMAT contains transformation matrix X^dag S X = 1
  ! NONOTMP is an array for storing AB in the product ABC...
  ! HORTHO contains X^dag H X

  INTEGER :: NONO_LWORK, NONO_LIWORK
  INTEGER :: NONZERO
  INTEGER, ALLOCATABLE  :: NONO_IWORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: NONO_WORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: UMAT(:,:), NONO_EVALS(:)
  REAL(LATTEPREC), ALLOCATABLE :: XMAT(:,:), SMAT(:,:), NONOTMP(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: ORTHOH(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: ORTHOHUP(:,:), ORTHOHDOWN(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HJJ(:), SH2(:,:)

  ! For the Pulay force = 2Tr[S^-1 H rho dS/dR ]

  REAL(LATTEPREC), ALLOCATABLE :: X2HRHO(:,:), SPINTMP(:,:)

END MODULE NONOARRAY
