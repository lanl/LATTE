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

SUBROUTINE READPPOTTAB

  USE CONSTANTS_MOD
  USE PPOTARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, MAXENTRY, NUMENTRY, N
  REAL(LATTEPREC) :: JUNK, P, QN, SIG, UN
  REAL(LATTEPREC), ALLOCATABLE :: U(:)
  CHARACTER(LEN=20) :: HD, HD1, HD2
  IF (EXISTERROR) RETURN

  OPEN(UNIT=14,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/ppots.dftb")

  READ(14,*) NOPPS

  ! Figure out array dimensions for allocation

  MAXENTRY = 0
  DO I = 1, NOPPS
     READ(14,*) HD1, HD2
     READ(14,*) NUMENTRY

     IF (NUMENTRY .GT. MAXENTRY) MAXENTRY = NUMENTRY

     DO J = 1, NUMENTRY
        READ(14,*) JUNK, JUNK
     ENDDO
  ENDDO

  REWIND(14)

  !  print*, MAXENTRY
  ! Now we can allocate

  ALLOCATE(PPELE1(NOPPS), PPELE2(NOPPS), PPR(MAXENTRY,NOPPS), &
       PPVAL(MAXENTRY,NOPPS))
  ALLOCATE(PPTABLENGTH(NOPPS), PPSPL(MAXENTRY,NOPPS))

  PPR = ZERO
  PPVAL = ZERO

  READ(14,*) NOPPS
  DO I = 1, NOPPS
     READ(14,*) PPELE1(I), PPELE2(I)
     READ(14,*) PPTABLENGTH(I)
     DO J = 1, PPTABLENGTH(I)
        READ(14,*) PPR(J,I), PPVAL(J,I)
     ENDDO
  ENDDO

  CLOSE(14)

  ! Let's get the cubic spline coefficients

  ALLOCATE(U(MAXENTRY))

  DO I = 1, NOPPS

     N = PPTABLENGTH(I)

     PPSPL(1,I) = ZERO
     U(1) = ZERO

     DO J = 2, N-1
        SIG = (PPR(J,I) - PPR(J-1,I))/(PPR(J+1,I) - PPR(J-1,I))
        P = SIG*PPSPL(J-1,I) + TWO
        PPSPL(J,I) = (SIG - ONE)/P
        U(J) = (SIX*((PPVAL(J+1,I) - PPVAL(J,I)) / &
             (PPR(J+1,I) - PPR(J,I)) - (PPVAL(J,I) - PPVAL(J-1,I)) &
             /(PPR(J,I) - PPR(J-1,I)))/(PPR(J+1,I)-PPR(J-1,I)) &
             - SIG*U(J-1))/P
     ENDDO

     QN = ZERO
     UN = ZERO

     PPSPL(N,I) = (UN - QN*U(N-1))/(QN*PPSPL(N-1, I) + ONE)

     DO K = N-1, 1, -1
        PPSPL(K,I) = PPSPL(K,I)*PPSPL(K+1,I) + U(K)
     ENDDO

  ENDDO

  DEALLOCATE(U)

  ! Test the interpolation

  !  CALL TABTEST

  RETURN

END SUBROUTINE READPPOTTAB
