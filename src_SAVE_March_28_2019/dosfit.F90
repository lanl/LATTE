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

SUBROUTINE DOSFIT

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, JUNK, ACC, II, COUNT
  INTEGER :: NSTEP, PICK, CPLOC(1), LOOPTARGET, AB
  REAL(LATTEPREC) :: EMIN, EMAX, MERIT, OLDMERIT, RN, ESTEP, EF
  REAL(LATTEPREC), ALLOCATABLE :: DFTDOS(:), DFTINTDOS(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBDOS(:), TBOLD(:), TBORIG(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBSCLOLD(:), TBSCLORIG(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBBEST(:)
  REAL(LATTEPREC), PARAMETER ::  ETA = 0.05, DOSWEIGHT = 1.0D0
  REAL(LATTEPREC) :: ENERGY, INTDOS, NUME, JUNKNUM, MINERR
  COMPLEX(LATTEPREC) :: CMPARG
  REAL(LATTEPREC), EXTERNAL :: GAUSSRN
  IF (EXISTERROR) RETURN

  MINERR = 1.0D6

  ! Read in the DOS from a VASP calculation

  OPEN(UNIT=20, STATUS="OLD", FILE="DOSCAR2fit.dat")

  READ(20,*) EMAX, EMIN, NSTEP, EF, JUNKNUM

  EMAX = EMAX - EF
  EMIN = EMIN - EF

  LOOPTARGET  = INT(BNDFIL*REAL(HDIM*NKTOT))

  ALLOCATE(DFTDOS(NSTEP), DFTINTDOS(NSTEP), TBDOS(NSTEP))
  ALLOCATE(TBOLD(INT2FIT), TBORIG(INT2FIT))
  ALLOCATE(TBSCLOLD(INT2FIT), TBSCLORIG(INT2FIT))
  ALLOCATE(TBBEST(INT2FIT))
  DO I = 1, NSTEP

     READ(20,*) JUNKNUM, DFTDOS(I), DFTINTDOS(I)

  ENDDO

  CLOSE(20)

  ! Initialize stuff


  IF (CONTROL .EQ. 1) THEN
     CALL ALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL ALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIALLOCATE
  ENDIF

  IF (ELECTRO .EQ. 1) THEN
     CALL ALLOCATECOULOMB
     CALL INITCOULOMB
  ENDIF

  !
  ! Allocate stuff for building the neighbor lists
  !

  CALL ALLOCATENEBARRAYS

  CALL NEBLISTS(0)  

  IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

  ! Compute DOS from LATTE and compare

  CALL INITRNG

  ACC = 0

  DO I = 1, INT2FIT
     TBORIG(I) = BOND(1,I)
     TBSCLORIG(I) = BOND(2,I)
  ENDDO

  NUME = ZERO
  DO I = 1, NATS
     NUME = NUME + ATOCC(ELEMPOINTER(I))
  ENDDO

  DO I = 1, INT2FIT
     TBBEST(I) = BOND(1,I)
  ENDDO



  DO II = 1, NFITSTEP


     ! Change set of TB parameters

     DO I = 1, INT2FIT

        TBOLD(I) = BOND(1,I)
        TBSCLOLD(I) = BOND(2,I)
     ENDDO

     ! Pick one to change


     IF (II .GT. 1) THEN

        CALL RANDOM_NUMBER(RN)

        PICK = INT(RN*REAL(INT2FIT)) + 1

        CALL RANDOM_NUMBER(RN)

        !        AB = INT(RN*2.0D0) + 1
        AB = 1
        CALL RANDOM_NUMBER(RN)

        ! Change by +/- 20%

        BOND(AB,PICK) = BOND(AB,PICK) * GAUSSRN(ONE, MCSIGMA)

        BOND(AB,PICK) = SIGN(BOND(AB,PICK), TBORIG(PICK))

        IF (ABS(BOND(1,PICK)) .LT. 0.01) &
             BOND(1,PICK) = SIGN(0.01D0, TBORIG(PICK))

        !        IF (BOND(1,PICK)/TBORIG(PICK) .LT. 0.5D0) THEN
        !           BOND(1,PICK) = 0.5D0*TBORIG(PICK)
        !        ELSEIF (BOND(1,PICK)/TBORIG(PICK) .GT. 1.5D0) THEN
        !           BOND(1,PICK) = 1.5D0*TBORIG(PICK)
        !        ENDIF

        !        IF (BOND(2,PICK)/TBSCLORIG(PICK) .LT. 0.5D0) THEN
        !           BOND(2,PICK) = 0.5D0*TBSCLORIG(PICK)
        !        ELSEIF (BOND(1,PICK)/TBORIG(PICK) .GT. 1.5D0) THEN
        !           BOND(2,PICK) = 1.5D0*TBSCLORIG(PICK)
        !        ENDIF


        CALL UNIVTAILCOEF(BOND(:,PICK))

     ENDIF

     ! Re-build cut-off tails

     ! Build H

     CALL KBLDNEWH

     IF (BASISTYPE .EQ. "NONORTHO") CALL KORTHOMYH

     IF (QFIT .EQ. 0) THEN

        CALL KDIAGMYH

        ! Find the chemical potential - we need the looptarget'th lowest 
        ! eigenvalue

        COUNT = 0
        DO I = 1, NKTOT
           DO J = 1, HDIM
              COUNT = COUNT + 1
              CPLIST(COUNT) = KEVALS(J,I)
           ENDDO
        ENDDO

        DO I = 1, LOOPTARGET
           CHEMPOT = MINVAL(CPLIST)
           CPLOC = MINLOC(CPLIST)
           CPLIST(CPLOC) = 1.0D12 ! We do this so we don't get this eigenvalue again on the next loop
        ENDDO

     ELSE


        !        IF (ELECTRO .EQ. 0) THEN
        !           CALL QNEUTRAL(0, 1)  ! Local charge neutrality
        !        ELSE
        !           CALL QCONSISTENCY(0,1) ! Self-consistent charges
        !        ENDIF


     ENDIF


     KEVALS = KEVALS - CHEMPOT

     ! Compute DOS

     TBDOS = ZERO

     ESTEP = (EMAX - EMIN)/REAL(NSTEP)

     DO I = 1, NSTEP

        ENERGY = EMIN + REAL(I-1)*ESTEP

        DO K = 1, NKTOT
           DO J = 1, HDIM

              CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))

              TBDOS(I) = TBDOS(I) - (ONE/PI)*AIMAG(CMPARG)

           ENDDO
        ENDDO

     ENDDO

     TBDOS = TBDOS/REAL(NKTOT)


     ! Normalize the DOS -> integral to Ef = Ne


     INTDOS = ZERO
     COUNT = 0
     DO I = 1, NSTEP

        ENERGY = EMIN + REAL(I-1)*ESTEP

        IF (ENERGY .LE. ZERO) THEN
           COUNT = COUNT + 1
           IF (MOD(I,2) .EQ. 0) THEN
              INTDOS = INTDOS + FOUR*TBDOS(I)
           ELSE
              INTDOS = INTDOS + TWO*TBDOS(I)
           ENDIF
        ENDIF

     ENDDO

     INTDOS = INTDOS - TBDOS(1) - TBDOS(COUNT)

     INTDOS = INTDOS*ESTEP/THREE

     TBDOS = TBDOS*NUME/INTDOS

     DO I = 1, NSTEP

        TBDOS(I) = ABS(TBDOS(I) - DFTDOS(I))

     ENDDO

     INTDOS = ZERO
     DO I = 1, NSTEP

        ENERGY = EMIN + REAL(I-1)*ESTEP

        IF (MOD(I,2) .EQ. 0) THEN
           INTDOS = INTDOS + FOUR*TBDOS(I)
        ELSE
           INTDOS = INTDOS + TWO*TBDOS(I)
        ENDIF

     ENDDO

     INTDOS = INTDOS - TBDOS(1) - TBDOS(NSTEP)

     INTDOS = INTDOS*ESTEP/THREE

     MERIT = INTDOS

     IF (II .EQ. 1) THEN
        OLDMERIT = MERIT
        MINERR = MERIT
     ENDIF


     IF (MERIT .LT. OLDMERIT) THEN

        ACC = ACC + 1
        OLDMERIT = MERIT

        IF (MERIT .LT. MINERR) THEN

           MINERR = MERIT
           DO I = 1, INT2FIT
              TBBEST(I) = BOND(1,I)
           ENDDO

        ENDIF

     ELSE

        CALL RANDOM_NUMBER(RN)

        IF ( EXP( -(MERIT - OLDMERIT)*MCBETA) .GT. RN ) THEN

           ACC = ACC + 1
           OLDMERIT = MERIT

        ELSE

           DO I = 1, INT2FIT

              BOND(1,I) = TBOLD(I)
              BOND(2,I) = TBSCLOLD(I)

           ENDDO


        ENDIF

     ENDIF

     WRITE(*,11) II, ACC, MERIT, OLDMERIT, &
          (BOND(1,J), J = 1, INT2FIT)
11   FORMAT(2I9, 2F12.6, 50F12.6) 

  ENDDO

  DO I = 1, INT2FIT
     BOND(1,I) = TBBEST(I)
  ENDDO


  DO I = 1, INT2FIT                                                 
     CALL UNIVTAILCOEF(BOND(:,I))                                   
  ENDDO

  CALL KBLDNEWH

  IF (BASISTYPE .EQ. "NONORTHO") CALL KORTHOMYH

  IF (QFIT .EQ. 0) THEN

     CALL KDIAGMYH

     ! Find the chemical potential - we need the looptarget'th lowest 
     ! eigenvalue

     COUNT = 0
     DO I = 1, NKTOT
        DO J = 1, HDIM
           COUNT = COUNT + 1
           CPLIST(COUNT) = KEVALS(J,I)
        ENDDO
     ENDDO

     DO I = 1, LOOPTARGET
        CHEMPOT = MINVAL(CPLIST)
        CPLOC = MINLOC(CPLIST)
        CPLIST(CPLOC) = 1.0D12 ! We do this so we don't get this eigenvalue again on the next loop
     ENDDO

  ELSE

     !        IF (ELECTRO .EQ. 0) THEN
     !           CALL QNEUTRAL(0, 1)  ! Local charge neutrality
     !        ELSE
     !           CALL QCONSISTENCY(0,1) ! Self-consistent charges
     !        ENDIF


  ENDIF

  KEVALS = KEVALS - CHEMPOT


  ! Build H                                                         


  ! Compute and normalize the DOS one more time

  TBDOS = ZERO

  ESTEP = (EMAX - EMIN)/REAL(NSTEP)

  DO I = 1, NSTEP

     ENERGY = EMIN + REAL(I-1)*ESTEP

     DO K = 1, NKTOT
        DO J = 1, HDIM

           CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))

           TBDOS(I) = TBDOS(I) - (ONE/PI)*AIMAG(CMPARG)

        ENDDO
     ENDDO

  ENDDO

  TBDOS = TBDOS/REAL(NKTOT)


  ! Normalize the DOS -> integral to Ef = Ne


  INTDOS = ZERO
  COUNT = 0
  DO I = 1, NSTEP

     ENERGY = EMIN + REAL(I-1)*ESTEP

     IF (ENERGY .LE. ZERO) THEN
        COUNT = COUNT + 1
        IF (MOD(I,2) .EQ. 0) THEN
           INTDOS = INTDOS + FOUR*TBDOS(I)
        ELSE
           INTDOS = INTDOS + TWO*TBDOS(I)
        ENDIF
     ENDIF

  ENDDO

  INTDOS = INTDOS - TBDOS(1) - TBDOS(COUNT)

  INTDOS = INTDOS*ESTEP/THREE

  TBDOS = TBDOS*NUME/INTDOS

  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="checkdosfit.dat")

  DO I = 1, NSTEP

     ENERGY = EMIN + REAL(I-1)*ESTEP

     WRITE(20,10) ENERGY, TBDOS(I), DFTDOS(I)

  ENDDO

  DO I = 1, INT2FIT

     PRINT*, "# ", BOND(1,I)
     WRITE(20,*) "# ", BOND(1,I)

  ENDDO

  CLOSE(20)

10 FORMAT(F12.3, 2G18.9)

  DEALLOCATE(TBDOS, DFTDOS, DFTINTDOS, TBSCLORIG, TBSCLOLD, TBOLD, TBORIG)
  DEALLOCATE(TBBEST)

  IF (CONTROL .EQ. 1) THEN
     !     CALL DEALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL DEALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIDEALLOCATE
  ENDIF

  IF (ELECTRO .EQ. 1)  CALL DEALLOCATECOULOMB

  IF (BASISTYPE .EQ. "NONORTHO") CALL DEALLOCATENONO

  CALL DEALLOCATENEBARRAYS


  RETURN

END SUBROUTINE DOSFIT
