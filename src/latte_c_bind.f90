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
!! \param NATS Number of atoms
!! \param COORDS Coordinates. Example: y-coordinate of atom 1 = COORDS(2,1)
!! \param TYPES An index for all the different atoms in the system. 
!! \param NTYPES Number of different elements in the system. 
!! \param MASSES Element masses for every different element of the system. 
!! \param XLO Lowest dimensions of the box.
!! \param XHI Highest dimensions of the box.
!! \param FORCES Forces for every atom as output.
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
SUBROUTINE LATTE_C_BIND (NATS,COORDS,TYPES,NTYPES,MASSES,XLO,XHI,FORCES) BIND (C, NAME="latte")

  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT
  USE LATTE_LIB    
   
  IMPLICIT NONE
  
  INTEGER(C_INT)                 :: NATS, NTYPES
  INTEGER(C_INT)                 ::  TYPES(NATS)
  REAL(C_DOUBLE)                 ::  COORDS(3,NTYPES), MASSES(NTYPES), XHI(3)
  REAL(C_DOUBLE)                 ::  XLO(3)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  FORCES(3, NATS)
   
  INTEGER                        ::  I
  
  CALL LATTE(NTYPES,TYPES,COORDS,MASSES,XLO,XHI,FORCES)
   
  RETURN
  
END SUBROUTINE LATTE_C_BIND       


