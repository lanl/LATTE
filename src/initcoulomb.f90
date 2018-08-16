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

SUBROUTINE INITCOULOMB

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: A2XA3(3), A3XA1(3), A1XA2(3)
  REAL(LATTEPREC) :: DOT, P, TIMERATIO, SQRTP
  IF (EXISTERROR) RETURN

  KECONST = 14.3996437701414*RELPERM
  TFACT  = 16.0/(5.0 * KECONST)

  IF (ELECMETH .EQ. 1) THEN

     !
     ! Get the coefficients for the cut-off tail for 1/R
     !

     CALL COULTAILCOEF

     COULCUT=ABS(COULCUT)
     COULCUT2 = COULCUT*COULCUT

  ELSE

     TWOPI = TWO*PI
     PI2 = PI*PI
     SQRTPI = SQRT(PI)
     EIGHTPI = EIGHT*PI

     ! First let's set up Ed's lattice_vecs

     LATTICEVECS = BOX

     ! Ed's bit:

     A2XA3(1) = LATTICEVECS(2,2)*LATTICEVECS(3,3) - &
          LATTICEVECS(2,3)*LATTICEVECS(3,2)
     A2XA3(2) = LATTICEVECS(2,3)*LATTICEVECS(3,1) - &
          LATTICEVECS(2,1)*LATTICEVECS(3,3)
     A2XA3(3) = LATTICEVECS(2,1)*LATTICEVECS(3,2) - &
          LATTICEVECS(2,2)*LATTICEVECS(3,1)

     DOT = LATTICEVECS(1,1)*A2XA3(1) + LATTICEVECS(1,2)*A2XA3(2) + &
          LATTICEVECS(1,3)*A2XA3(3)

     RECIPVECS(1,1) = TWOPI*A2XA3(1)/DOT
     RECIPVECS(1,2) = TWOPI*A2XA3(2)/DOT
     RECIPVECS(1,3) = TWOPI*A2XA3(3)/DOT

     A3XA1(1) = LATTICEVECS(3,2)*LATTICEVECS(1,3) - &
          LATTICEVECS(3,3)*LATTICEVECS(1,2)
     A3XA1(2) = LATTICEVECS(3,3)*LATTICEVECS(1,1) - &
          LATTICEVECS(3,1)*LATTICEVECS(1,3)
     A3XA1(3) = LATTICEVECS(3,1)*LATTICEVECS(1,2) - &
          LATTICEVECS(3,2)*LATTICEVECS(1,1)

     DOT = LATTICEVECS(2,1)*A3XA1(1) + LATTICEVECS(2,2)*A3XA1(2) + &
          LATTICEVECS(2,3)*A3XA1(3)

     RECIPVECS(2,1) = TWOPI*A3XA1(1)/DOT
     RECIPVECS(2,2) = TWOPI*A3XA1(2)/DOT
     RECIPVECS(2,3) = TWOPI*A3XA1(3)/DOT

     A1XA2(1) = LATTICEVECS(1,2)*LATTICEVECS(2,3) - &
          LATTICEVECS(1,3)*LATTICEVECS(2,2)
     A1XA2(2) = LATTICEVECS(1,3)*LATTICEVECS(2,1) - &
          LATTICEVECS(1,1)*LATTICEVECS(2,3)
     A1XA2(3) = LATTICEVECS(1,1)*LATTICEVECS(2,2) - &
          LATTICEVECS(1,2)*LATTICEVECS(2,1)

     DOT = LATTICEVECS(3,1)*A1XA2(1) + LATTICEVECS(3,2)*A1XA2(2) + &
          LATTICEVECS(3,3)*A1XA2(3)

     RECIPVECS(3,1) = TWOPI*A1XA2(1)/DOT
     RECIPVECS(3,2) = TWOPI*A1XA2(2)/DOT
     RECIPVECS(3,3) = TWOPI*A1XA2(3)/DOT

     ! Calculate the cell volume

     COULVOL = DOT

     P = -LOG(COULACC)
     SQRTP = SQRT(P)

     IF (COULCUT .GT. ZERO) THEN

        CALPHA = SQRTP/COULCUT
        COULCUT2 = COULCUT*COULCUT
        KCUTOFF = TWO*CALPHA*SQRTP
        KCUTOFF2 = KCUTOFF*KCUTOFF
        CALPHA2 = CALPHA*CALPHA
        FOURCALPHA2 = FOUR*CALPHA2

     ELSE

        !
        ! Automatically determining the optimal real space
        ! cut-off if on input COUTCUT < 0. This is Sanville's code
        !

        !   TIMERATIO = 50.0
        !   TIMERATIO = 1.00D0
        TIMERATIO = 3.50D0 

        CALPHA = SQRTPI*((TIMERATIO * REAL(NATS) / (COULVOL*COULVOL))**(ONE/SIX))
        COULCUT = SQRTP/CALPHA

        !        PRINT*, "COULCUT =", COULCUT
        IF (COULCUT .GT. FIVE*TEN) THEN

           COULCUT = FIVE*TEN
           CALPHA = SQRTP/COULCUT

        ENDIF

        COULCUT2 = COULCUT*COULCUT
        KCUTOFF = TWO*CALPHA*SQRTP
        KCUTOFF2 = KCUTOFF*KCUTOFF
        CALPHA2 = CALPHA*CALPHA
        FOURCALPHA2 = FOUR*CALPHA2

        ! Taking this bit from Coulomb Ewald so we don't have to
        ! recompute every time:

     ENDIF

     LMAX = INT(KCUTOFF / SQRT(RECIPVECS(1,1)*RECIPVECS(1,1) + &
          RECIPVECS(1,2)*RECIPVECS(1,2) + RECIPVECS(1,3)*RECIPVECS(1,3)))

     MMAX = INT(KCUTOFF / SQRT(RECIPVECS(2,1)*RECIPVECS(2,1) + &
          RECIPVECS(2,2)*RECIPVECS(2,2) + RECIPVECS(2,3)*RECIPVECS(2,3)))

     NMAX = INT(KCUTOFF / SQRT(RECIPVECS(3,1)*RECIPVECS(3,1) + &
          RECIPVECS(3,2)*RECIPVECS(3,2) + RECIPVECS(3,3)*RECIPVECS(3,3)))


  ENDIF

  CALL GET_K_LISTS(RECIPVECS)

END SUBROUTINE INITCOULOMB
