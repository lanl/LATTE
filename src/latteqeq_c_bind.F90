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

!> Subroutine for Latte-C++ interfacing.

SUBROUTINE LATTEQEQ_C_BIND (FLAGS, NATS, COORDS, TYPES, NTYPES, MASSES, XLO, &
     XHI, XY, XZ, YZ, FORCES, MAXITER, VENERG, &
     VEL, DT, VIRIAL_INOUT, CURRENTSTEP, &!QXLBO, KK, 
     GRADX, &
     NEWSYSTEM, EXISTERROR) BIND (C, NAME="latteqeq")

  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_BOOL
  USE qeq_LIB

  IMPLICIT NONE

  INTEGER(C_INT)                 ::  CURRENTSTEP, NATS, NTYPES, MAXITER
  INTEGER(C_INT)                 ::  TYPES(NATS), FLAGS(5)
  !REAL(C_DOUBLE)                 ::  QXLBO(NATS,8), KK(NATS, NATS)
  REAL(C_DOUBLE)                 ::  GRADX(NATS*3+1, NATS)
  REAL(C_DOUBLE)                 ::  COORDS(3,NATS), MASSES(NTYPES), XHI(3)
  REAL(C_DOUBLE)                 ::  XLO(3), EKIN, VENERG, DT
  REAL(C_DOUBLE)                 ::  XY, XZ, YZ
  REAL(C_DOUBLE), INTENT(INOUT)  ::  FORCES(3, NATS), VEL(3, NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  VIRIAL_INOUT(6)
  LOGICAL(C_BOOL)                ::  EXISTERROR
  INTEGER(C_INT), INTENT(INOUT)  ::  NEWSYSTEM

  CALL LATTEQEQ(NTYPES, TYPES, COORDS, MASSES, XLO, XHI, XY, XZ, YZ, FORCES, &
       MAXITER, VENERG, VEL, DT, VIRIAL_INOUT, CURRENTSTEP, & !QXLBO, KK, &
       GRADX, NEWSYSTEM, EXISTERROR)

  RETURN

END SUBROUTINE LATTEQEQ_C_BIND

!> Function for Latte-C++ interfacing.
!! \return ABIVERSION integer representing the date of the last change
!!          to the C/C++ interface (e.g. 20180221)
!!
!! \brief This function will be used prior to calling the LATTE library
!!        to allow the calling code to ensure the linked library version is compatible.
!!
INTEGER(C_INT) FUNCTION LATTEQEQ_C_ABIVERSION()  BIND (C, NAME="latteqeq_abiversion")
  USE ISO_C_BINDING, ONLY: C_INT
  USE qeq_LIB,     ONLY: LATTEQEQ_ABIVERSION
  IMPLICIT NONE

  LATTEQEQ_C_ABIVERSION = LATTEQEQ_ABIVERSION
  RETURN

END FUNCTION LATTEQEQ_C_ABIVERSION
