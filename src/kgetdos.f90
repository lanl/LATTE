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

SUBROUTINE KGETDOS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE MYPRECISION
  USE KSPACEARRAY
  USE MDARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, COUNT
  REAL(LATTEPREC), PARAMETER :: EWINDOW = 0.05, ETA = 0.05
  REAL(LATTEPREC) :: INTDOS, ENERGY, NUME, DOSMIN, DOSMAX
  INTEGER :: NDOSBINS, BINID, NUMORB
  REAL(LATTEPREC), ALLOCATABLE :: DOS(:)
  COMPLEX(LATTEPREC) :: CMPARG
  IF (EXISTERROR) RETURN

  DOSMAX = MAXVAL(KEVALS) + 10.0*EWINDOW
  DOSMIN = MINVAL(KEVALS) - 10.0*EWINDOW

  NDOSBINS = INT((DOSMAX - DOSMIN)/EWINDOW)

  ALLOCATE(DOS(NDOSBINS))

  DOS = ZERO

  DO I = 1, NDOSBINS

     ENERGY = DOSMIN + REAL(I-1)*EWINDOW

     DO K = 1, NKTOT
        DO J = 1, HDIM

           CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))

           DOS(I) = DOS(I) - (ONE/PI)*AIMAG(CMPARG)

        ENDDO
     ENDDO

  ENDDO

  DOS = DOS/REAL(NKTOT)

  ! To normalize, lets integrate the dos

  INTDOS = ZERO
  COUNT = 0
  DO I = 1, NDOSBINS

     ENERGY = DOSMIN + REAL(I-1)*EWINDOW

     IF (ENERGY .LE. CHEMPOT) THEN
        COUNT = COUNT + 1
        IF (MOD(I,2) .EQ. 0) THEN
           INTDOS = INTDOS + FOUR*DOS(I)
        ELSE
           INTDOS = INTDOS + TWO*DOS(I)
        ENDIF
     ENDIF

  ENDDO

  INTDOS = INTDOS - DOS(1) - DOS(COUNT)

  INTDOS = INTDOS*EWINDOW/THREE

  ! Let's figure out what the integrated DOS should be
  ! In VASP the integral of the DOS up to the chemical 
  ! potential = the total number of electrons. Let's do the 
  ! same.

  NUME = ZERO
  DO I = 1, NATS
     NUME = NUME + ATOCC(ELEMPOINTER(I))
  ENDDO

  DOS = DOS*NUME/INTDOS


  OPEN(UNIT=50, STATUS="UNKNOWN", FILE="DoS.dat")

  DO I = 1, NDOSBINS

     WRITE(50, 10) DOSMIN + REAL(I-1)*EWINDOW - CHEMPOT, DOS(I)

  ENDDO

10 FORMAT(2F12.6)

  CLOSE(50)

  DEALLOCATE(DOS)

  RETURN

END SUBROUTINE KGETDOS

