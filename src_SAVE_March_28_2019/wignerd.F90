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

FUNCTION WIGNERD(L, M, MP, COSBETA)

  ! Builds Wigner d function 
  ! notation conforms to that in PRB 72 165107 (2005), eq. (9)

  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: BETA, PREF, SUM, NUMER, DENOM, WIGNERD
  REAL(LATTEPREC) :: COSBETA
  INTEGER, EXTERNAL :: FACTORIAL
  INTEGER :: K, L, M, MP

  !  PREF = TWO ** (-L) * REAL( (-1)**(L - MP)) * &
  !       SQRT( REAL(FACTORIAL(L + M) * FACTORIAL(L - M) * &
  !       FACTORIAL(L + MP) * FACTORIAL(L - MP)) )

  PREF =  REAL( (-1)**(L - MP)) * &
       SQRT( REAL(FACTORIAL(L + M) * FACTORIAL(L - M) * &
       FACTORIAL(L + MP) * FACTORIAL(L - MP)) ) / REAL(2**L)

  SUM = ZERO

  DO K = MAX(0, - M - MP), MIN(L - M, L - MP)

     NUMER = REAL( (-1)**K ) * &
          (ONE - COSBETA) ** (L - K - HALF * (M + MP)) * &
          (ONE + COSBETA) ** (K + HALF * (M + MP))

     DENOM = REAL(FACTORIAL(K) * FACTORIAL(L - M - K) * &
          FACTORIAL(L - MP - K) * FACTORIAL(M + MP + K))

     SUM = SUM + NUMER / DENOM

  ENDDO

  WIGNERD = PREF * SUM

  RETURN 

END FUNCTION WIGNERD
