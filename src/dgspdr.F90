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

FUNCTION DGSPDR(I, J, L1, L2, M, MAGR)

  ! Derivative of Goodwin-Skinner-Pettifor function
  ! Note that -d/dR is returned

  USE GSPARRAY
  USE MYPRECISION
  USE SETUPARRAY

  IMPLICIT NONE

  INTEGER :: IC, IGS, I, J, L1, L2, M, IP1, IP2
  REAL(LATTEPREC) :: MAGR, DGSPDR, DC
  REAL(LATTEPREC) :: RMINUSR1, RMINUSR12, RMINUSR14
  REAL(LATTEPREC), EXTERNAL :: GSP
  CHARACTER (LEN = 3) :: IGLTYPE
  LOGICAL BINGO

  ! can't test directly on L values because basis strings always list
  ! lower L values first
  IF (L1 > L2) THEN
     IP1 = L2
     IP2 = L1
  ELSE
     IP1 = L1
     IP2 = L2
  ENDIF

  ! build basis string from L and M values - pure hackery
  SELECT CASE(IP1)
  CASE(0)
     IGLTYPE = "s"
  CASE(1)
     IGLTYPE = "p"
  CASE(2)
     IGLTYPE = "d"
  CASE(3)
     IGLTYPE = "f"
  END SELECT
  SELECT CASE(IP2)
  CASE(0)
     IGLTYPE = TRIM(IGLTYPE)//"s"
  CASE(1)
     IGLTYPE = TRIM(IGLTYPE)//"p"
  CASE(2)
     IGLTYPE = TRIM(IGLTYPE)//"d"
  CASE(3)
     IGLTYPE = TRIM(IGLTYPE)//"f"
  END SELECT
  SELECT CASE(M)
  CASE(0)
     IGLTYPE = TRIM(IGLTYPE)//"s"
  CASE(1)
     IGLTYPE = TRIM(IGLTYPE)//"p"
  CASE(2)
     IGLTYPE = TRIM(IGLTYPE)//"d"
  CASE(3)
     IGLTYPE = TRIM(IGLTYPE)//"f"
  END SELECT

  IC = 0
  BINGO = .FALSE.
  DO WHILE (.NOT. BINGO)
     IC = IC + 1
     IF (((ATELE(I) .EQ. ELE1(IC) .AND. ATELE(J) .EQ. ELE2(IC)) .OR. &
          (ATELE(J) .EQ. ELE1(IC) .AND. ATELE(I) .EQ. ELE2(IC))) .AND. &
          (IGLTYPE .EQ. BTYPE(IC))) THEN
        IGS = IC
        BINGO = .TRUE.
     ENDIF
  ENDDO

  IF (MAGR .LE. GSPR1(IGS)) THEN

     DGSPDR = GSPN(IGS) * GSP(I, J, L1, L2, M, MAGR) * (GSPNC(IGS)* &
          ((MAGR/GSPRC(IGS))**GSPNC(IGS)) + ONE) / MAGR

  ELSEIF (MAGR .LT. GSPRCUT(IGS)) THEN

     RMINUSR1 = MAGR - GSPR1(IGS)
     RMINUSR12 = RMINUSR1*RMINUSR1
     RMINUSR14 = RMINUSR12*RMINUSR12

     DGSPDR = MINUSONE * (B(2,IGS) + TWO * B(3,IGS) * RMINUSR1 + &
          THREE * B(4,IGS) * RMINUSR12 + FOUR * B(5,IGS) * RMINUSR12 * &
          RMINUSR1 + FIVE * B(6,IGS) * RMINUSR14)
     DGSPDR = DGSPDR * HR0(IGS)

     ! note that permutation symmetry when GSPR1 < MAGR < GSPRCUT 
     ! is accounted for in GSP
     IF (L1 > L2) THEN
        IF (MOD(L1 + L2, 2) /= 0) DGSPDR = MINUSONE * DGSPDR
     ENDIF

  ELSE

     DGSPDR = ZERO

  ENDIF

  RETURN

END FUNCTION DGSPDR
