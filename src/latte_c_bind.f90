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
!! \param FLAGS Different control flags that can be passed to LATTE (not in use yet)
!! \param NATS Number of atoms
!! \param COORDS Coordinates. Example: y-coordinate of atom 1 = COORDS(2,1)
!! \param TYPES An index for all the different atoms in the system.
!! \param NTYPES Number of different elements in the system
!! \param MASSES Element masses for every different element of the system.
!! \param XLO Lowest dimensions of the box
!! \param XHI Highest dimensions of the box
!! \param XY Tilt factor.
!! \param XZ Tilt factor.
!! \param YZ Tilt factor. The lattice vectors are constructed as:
!! a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo).
!! \param FORCES Forces for every atom as output.
!! \param MAXITER Latte MAXITER keyword. If MAXITER = -1, only the Forces are computed.
!!        If MAXITER = 0, MAXITER is read from latte.in file.
!!        IF MAXITER > 0, MAXITER is passed trough the library call.
!! \param VENERG This is the potential Energy that is given back from latte to the hosting code.
!! \param VEL Velocities passed to latte.
!! \param DT integration step passed to latte.
!! \param DT integration step passed to latte.
!! \param VIRIAL_INOUT Components of the second virial coefficient
!! \param NEWSYSTEM Tells LATTE if a new system is passed.
!! \param EXISTERROR Returns an error flag (.true.) to the hosting code.
!!
!! \brief This routine will be used load call latte_lib from a C/C++ program:
!!
!! \brief Note: To get the mass of atom 3 we do:
!! \verbatim
!!      MASS(TYPES(3))
!! \endverbatim
!!
!! \brief Note: To get the lattice vectors as formated in LATTE we do:
!! \verbatim
!!      BOX(1,1) = XHI(1) - XLO(1); ...
!! \endverbatim
!!
!! \brief Note: All units are LATTE units by default.
!! See https://github.com/losalamos/LATTE/blob/master/Manual/LATTE_manual.pdf
!!
SUBROUTINE LATTE_C_BIND (FLAGS, NATS, COORDS, TYPES, NTYPES, MASSES, XLO &
     , XHI, XY, XZ, YZ, FORCES, MAXITER, VENERG, &
     VEL, DT, VIRIAL_INOUT, NEWSYSTEM, EXISTERROR) BIND (C, NAME="latte")

  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_BOOL
  USE LATTE_LIB

  IMPLICIT NONE

  INTEGER(C_INT)                 ::  NATS, NTYPES, MAXITER
  INTEGER(C_INT)                 ::  TYPES(NATS), FLAGS(5)
  REAL(C_DOUBLE)                 ::  COORDS(3,NATS), MASSES(NTYPES), XHI(3)
  REAL(C_DOUBLE)                 ::  XLO(3), EKIN, VENERG, DT
  REAL(C_DOUBLE)                 ::  XY, XZ, YZ
  REAL(C_DOUBLE), INTENT(INOUT)  ::  FORCES(3, NATS), VEL(3, NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  VIRIAL_INOUT(6)
  LOGICAL(C_BOOL)                ::  EXISTERROR
  INTEGER(C_INT), INTENT(INOUT)  ::  NEWSYSTEM

  CALL LATTE(NTYPES, TYPES, COORDS, MASSES, XLO, XHI, XY, XZ, YZ, FORCES, &
       MAXITER, VENERG, VEL, DT, VIRIAL_INOUT, NEWSYSTEM, EXISTERROR)

  RETURN

END SUBROUTINE LATTE_C_BIND

!> Function for Latte-C++ interfacing.
!! \return ABIVERSION integer representing the date of the last change
!!          to the C/C++ interface (e.g. 20180221)
!!
!! \brief This function will be used prior to calling the LATTE library
!!        to allow the calling code to ensure the linked library version is compatible.
!!
INTEGER(C_INT) FUNCTION LATTE_C_ABIVERSION()  BIND (C, NAME="latte_abiversion")
  USE ISO_C_BINDING, ONLY: C_INT
  USE LATTE_LIB,     ONLY: LATTE_ABIVERSION
  IMPLICIT NONE

  LATTE_C_ABIVERSION = LATTE_ABIVERSION
  RETURN

END FUNCTION LATTE_C_ABIVERSION
