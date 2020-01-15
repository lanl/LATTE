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

  USE WIGNERARRAY
  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: BETA, PREF, NUMER, DENOM, WIGNERD
  REAL(LATTEPREC) :: COSBETA, POWER1, POWER2, POWTMP
  INTEGER :: K, L, M, MP


  IF (ABS(MP) .GT. L) THEN

!     PRINT*, L+M, L-M, L+MP, L-MP
     WIGNERD = ZERO

  ELSE

     PREF = (REAL(MINUSONEPOW( L - MP ))/REAL(TWOPOW(L))) * &
          SQRTFACT(L+M)*SQRTFACT(L-M)*SQRTFACT(L+MP)*SQRTFACT(L-MP)

     WIGNERD = ZERO
  
     DO K = MAX(0, - M - MP), MIN(L - M, L - MP)

        POWTMP = REAL(M + MP)/TWO
        POWER1 = REAL(L - K) - POWTMP
        POWER2 = REAL(K) + POWTMP

!        PRINT*, POWER1, POWER2

        WIGNERD = WIGNERD + REAL( MINUSONEPOW(K) ) * &
             ((ONE - COSBETA) ** POWER1) * ((ONE + COSBETA) ** POWER2) / &
             REAL(FACTORIAL(K) * FACTORIAL(L - M - K) * &
             FACTORIAL(L - MP - K) * FACTORIAL(M + MP + K))
     ENDDO

     WIGNERD = PREF * WIGNERD

  ENDIF

  RETURN 

END FUNCTION WIGNERD
