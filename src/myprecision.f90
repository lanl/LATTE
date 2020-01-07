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

MODULE MYPRECISION

  IMPLICIT NONE
  SAVE

  !
  ! The precision for the calculation (double = 8, real = 4)
  ! Selected during the preprocessing of myprecision.F
  !

#ifdef DOUBLEPREC	
  INTEGER, PARAMETER :: LATTEPREC = 8
#elif defined(SINGLEPREC)
  INTEGER, PARAMETER :: LATTEPREC = 4
#endif

  !
  ! In the case of single precision, we assume (uh uh!) that the compiler
  ! will round our doubles down to singles for us...
  !

  REAL(LATTEPREC), PARAMETER :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0
  REAL(LATTEPREC), PARAMETER :: THREE = 3.0D0, FOUR = 4.0D0, FIVE = 5.0D0
  REAL(LATTEPREC), PARAMETER :: SIX = 6.0D0, SEVEN = 7.0D0, EIGHT = 8.0D0
  REAL(LATTEPREC), PARAMETER :: NINE = 9.0D0, TEN = 10.0D0, ELEVEN = 11.0D0
  REAL(LATTEPREC), PARAMETER :: TWELVE = 12.0D0, FOURTEEN = 14.0D0
  REAL(LATTEPREC), PARAMETER :: FIFTEEN = 15.0D0
  REAL(LATTEPREC), PARAMETER :: SIXTEEN = 16.0D0, TWENTY = 20.0D0
  REAL(LATTEPREC), PARAMETER :: TWENTYFOUR = 24.0D0, TWENTYSIX = 26.0D0
  REAL(LATTEPREC), PARAMETER :: FORTYEIGHT = 48.0D0
  REAL(LATTEPREC), PARAMETER :: HALF = 0.5D0, THIRD = 1.0D0/3.0D0
  REAL(LATTEPREC), PARAMETER :: QUARTER = 0.25D0, THREEQUART = 0.75D0
  REAL(LATTEPREC), PARAMETER :: MINUSONE = -1.0D0, THOUSAND = 1000.0D0
  REAL(LATTEPREC), PARAMETER :: SQRT2 = SQRT(2.0D0)
  COMPLEX(LATTEPREC), PARAMETER :: RE = (ONE, ZERO), IM = (ZERO,ONE)

  ! Short successions lists
  INTEGER :: FACTORIALLIST(0:10)
  REAL(LATTEPREC) :: MINUSONEPOWERS(0:10)

END MODULE MYPRECISION
