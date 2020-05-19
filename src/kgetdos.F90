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
  REAL(LATTEPREC), PARAMETER :: EWINDOW = 0.05D0, ETA = 0.05D0
  REAL(LATTEPREC) :: INTDOS, ENERGY, NUME, DOSMIN, DOSMAX
  REAL(LATTEPREC) :: INTDOSUP, INTDOSDOWN
  INTEGER :: NDOSBINS, BINID, NUMORB
  REAL(LATTEPREC), ALLOCATABLE :: DOS(:), DOSUP(:), DOSDOWN(:)
  COMPLEX(LATTEPREC) :: CMPARG
  IF (EXISTERROR) RETURN

  IF (SPINON .EQ. 0) THEN
     DOSMAX = MAXVAL(KEVALS) + 10.0D0*EWINDOW
     DOSMIN = MINVAL(KEVALS) - 10.0D0*EWINDOW
  ELSE
     DOSMAX = MAX(MAXVAL(KEVALSUP), MAXVAL(KEVALSDOWN)) + 10.0D0*EWINDOW
     DOSMIN = MIN(MINVAL(KEVALSUP), MINVAL(KEVALSDOWN)) - 10.0D0*EWINDOW
  ENDIF

  NDOSBINS = INT((DOSMAX - DOSMIN)/EWINDOW)

  IF (SPINON .EQ. 0) THEN
     ALLOCATE(DOS(NDOSBINS))
     DOS = ZERO
  ELSE
     ALLOCATE(DOSUP(NDOSBINS), DOSDOWN(NDOSBINS))
     DOSUP = ZERO
     DOSDOWN = ZERO
  ENDIF


  DO I = 1, NDOSBINS

     ENERGY = DOSMIN + REAL(I-1)*EWINDOW

     IF (SPINON .EQ. 0) THEN

        DO K = 1, NKTOT
           DO J = 1, HDIM

              CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))

              DOS(I) = DOS(I) - (ONE/PI)*AIMAG(CMPARG)

           ENDDO
        ENDDO
        
        DOS(I) = DOS(I)/REAL(NKTOT)

     ELSE

        DO K = 1, NKTOT
           DO J = 1, HDIM

              CMPARG = ONE/(ENERGY - KEVALSUP(J,K) + CMPLX(ZERO,ETA))

              DOSUP(I) = DOSUP(I) - (ONE/PI)*AIMAG(CMPARG)

              CMPARG = ONE/(ENERGY - KEVALSDOWN(J,K) + CMPLX(ZERO,ETA))

              DOSDOWN(I) = DOSDOWN(I) - (ONE/PI)*AIMAG(CMPARG)

           ENDDO
        ENDDO
        
        DOSUP(I) = DOSUP(I)/REAL(NKTOT)
        DOSDOWN(I) = DOSDOWN(I)/REAL(NKTOT)

     ENDIF

  ENDDO


  ! To normalize, lets integrate the dos

  NUME = ZERO
  DO I = 1, NATS
     NUME = NUME + ATOCC(ELEMPOINTER(I))
  ENDDO


  IF (SPINON .EQ. 0) THEN
     
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

     DOS = DOS*NUME/INTDOS

  ELSE

     INTDOSUP = ZERO
     INTDOSDOWN = ZERO

     COUNT = 0
     DO I = 1, NDOSBINS

        ENERGY = DOSMIN + REAL(I-1)*EWINDOW

        IF (ENERGY .LE. CHEMPOT) THEN
           COUNT = COUNT + 1
           IF (MOD(I,2) .EQ. 0) THEN
              INTDOS = INTDOS + FOUR*(DOSUP(I) + DOSDOWN(I))
!              INTDOSUP = INTDOSUP + FOUR*DOSUP(I)
!              INTDOSDOWN = INTDOSDOWN + FOUR*DOSDOWN(I)
           ELSE
              INTDOS = INTDOS + TWO*(DOSUP(I) + DOSDOWN(I))
!              INTDOSUP = INTDOSUP + TWO*DOSUP(I)
!              INTDOSDOWN = INTDOSDOWN + TWO*DOSDOWN(I)
           ENDIF
        ENDIF
        
     ENDDO

     INTDOS = INTDOS - (DOSUP(1) + DOSDOWN(1) + DOSUP(COUNT)+DOSDOWN(COUNT))
     
!     INTDOSUP = INTDOSUP - DOSUP(1) - DOSUP(COUNT)

!     INTDOSUP = INTDOSUP*EWINDOW/THREE

!     INTDOSDOWN = INTDOSDOWN - DOSDOWN(1) - DOSDOWN(COUNT)

!     INTDOSDOWN = INTDOSDOWN*EWINDOW/THREE

!     DOSUP = DOSUP*NUME/(INTDOSUP*TWO)
!     DOSDOWN = DOSDOWN*NUME/(INTDOSDOWN*TWO)

     INTDOS = INTDOS*EWINDOW/THREE

     DOSUP = DOSUP*NUME/(INTDOS*TWO)
     DOSDOWN = DOSDOWN*NUME/(INTDOS*TWO)

     
  ENDIF

  OPEN(UNIT=50, STATUS="UNKNOWN", FILE="DoS.dat")

  DO I = 1, NDOSBINS

     IF (SPINON .EQ. 0) THEN
        WRITE(50, 10) DOSMIN + REAL(I-1)*EWINDOW - CHEMPOT, DOS(I)
     ELSE
        WRITE(50, 10) DOSMIN + REAL(I-1)*EWINDOW - CHEMPOT, DOSUP(I), DOSDOWN(I)
     ENDIF

  ENDDO

10 FORMAT(3F12.6)

  CLOSE(50)

  IF (SPINON .EQ. 0) THEN
     DEALLOCATE(DOS)
  ELSE 
     DEALLOCATE(DOSUP, DOSDOWN)
  ENDIF

  RETURN

END SUBROUTINE KGETDOS

