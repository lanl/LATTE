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

SUBROUTINE MOFIT

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
  REAL(LATTEPREC), ALLOCATABLE :: DFTMO(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBDOS(:), TBOLD(:), TBORIG(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBSCLOLD(:), TBSCLORIG(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBBEST(:), TBBESTSCL(:)
  REAL(LATTEPREC) :: ENERGY, INTDOS, NUME, JUNKNUM, MINERR
  REAL(LATTEPREC), EXTERNAL :: GAUSSRN
  IF (EXISTERROR) RETURN

  ! Read in the DOS from a VASP calculation

  OPEN(UNIT=20, STATUS="OLD", FILE="MO2fit.dat")

  ALLOCATE(DFTMO(HDIM))
  ALLOCATE(TBOLD(INT2FIT), TBORIG(INT2FIT))
  ALLOCATE(TBSCLOLD(INT2FIT), TBSCLORIG(INT2FIT))
  ALLOCATE(TBBEST(INT2FIT), TBBESTSCL(INT2FIT))

  DO I = 1, HDIM

     READ(20,*) J, DFTMO(I), JUNKNUM

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
     TBBEST(I) = BOND(1,I)
     TBBESTSCL(I) = BOND(2,I)
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

        BOND(AB,PICK) = BOND(AB,PICK) * GAUSSRN(1.0, MCSIGMA)


        IF (ABS(BOND(1,PICK)) .LT. 0.01) &
             BOND(1,PICK) = SIGN(0.01D0, TBORIG(PICK))

        !        IF (BOND(1,PICK)/TBORIG(PICK) .LT. 0.5D0) THEN
        !           BOND(1,PICK) = 0.5D0*TBORIG(PICK)
        !        ELSEIF (BOND(1,PICK)/TBORIG(PICK) .GT. 1.5D0) THEN
        !           BOND(1,PICK) = 1.5D0*TBORIG(PICK)
        !        ENDIF

        !        IF (BOND(1,4) .GT. BOND(1,3)) BOND(1,3) = BOND(1,4)

        !        IF (ABS(BOND(1,3)/BOND(1,4)) .LT. 2.0D0) THEN
        !           BOND(1,3) = 2.0D0*BOND(1,4)
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

     CALL BLDNEWHS_SP

     IF (QFIT .EQ. 0) THEN

        CALL DIAGMYH

     ELSE

        IF (ELECTRO .EQ. 0) THEN
           CALL QNEUTRAL(0, 1)  ! Local charge neutrality
        ELSE
           CALL QCONSISTENCY(0,1) ! Self-consistent charges
        ENDIF

     ENDIF

     ! Compute errors

     MERIT = ZERO

     DO I = 1, HDIM

        MERIT = MERIT + (EVALS(I) - DFTMO(I))*(EVALS(I) - DFTMO(I))

     ENDDO

     IF (II .EQ. 1) OLDMERIT = MERIT
     IF (II .EQ. 1) MINERR = MERIT

     IF (MERIT .LT. OLDMERIT) THEN

        ACC = ACC + 1
        OLDMERIT = MERIT

        WRITE(*,11) II, ACC, MERIT, OLDMERIT, &
             (BOND(1,J), BOND(2,J), J = 1, INT2FIT)
11      FORMAT(2I9, 2F12.6, 50F12.6)         

        IF (MERIT .LT. MINERR) THEN
           MINERR = MERIT
           DO I = 1, INT2FIT
              TBBEST(I) = BOND(1,I)
              TBBESTSCL(I) = BOND(2,I)
           ENDDO
        ENDIF


     ELSE

        CALL RANDOM_NUMBER(RN)

        IF ( EXP( -(MERIT - OLDMERIT)*MCBETA) .GT. RN ) THEN

           ACC = ACC + 1
           OLDMERIT = MERIT

           WRITE(*,11) II, ACC, MERIT, OLDMERIT, &
                (BOND(1,J), BOND(2,J), J = 1, INT2FIT)

        ELSE

           DO I = 1, INT2FIT

              BOND(1,I) = TBOLD(I)
              BOND(2,I) = TBSCLOLD(I)

           ENDDO


        ENDIF

     ENDIF

     !     WRITE(*,11) II, ACC, MERIT, OLDMERIT, &
     !          (BOND(1,J), BOND(2,J), J = 1, INT2FIT)
     !     11 FORMAT(2I9, 2F12.6, 50F12.6) 

  ENDDO

  PRINT*, "MINERR = ", MINERR

  DO I = 1, INT2FIT
     BOND(1,I) = TBBEST(I)
     BOND(2,I) = TBBESTSCL(I)
  ENDDO


  DO I = 1, INT2FIT                                                 
     CALL UNIVTAILCOEF(BOND(:,I))                                   
  ENDDO


  CALL BLDNEWHS_SP

  IF (QFIT .EQ. 0) THEN

     CALL DIAGMYH

  ELSE

     IF (ELECTRO .EQ. 0) THEN
        CALL QNEUTRAL(0, 1)  ! Local charge neutrality
     ELSE
        CALL QCONSISTENCY(0,1) ! Self-consistent charges
     ENDIF


  ENDIF


  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="checkmofit.dat")

  DO I = 1, HDIM

     WRITE(20,10) EVALS(I), DFTMO(I)

  ENDDO

  DO I = 1, INT2FIT

     PRINT*, "# ", BOND(1,I)
     WRITE(20,*) "# ", BOND(1,I)

  ENDDO

  WRITE(20,20) (BOND(1,I), I = 1, INT2FIT)

  CLOSE(20)

10 FORMAT(F12.3, 2G18.9)
20 FORMAT(100G18.9)

  DEALLOCATE(DFTMO, TBSCLORIG, TBSCLOLD, TBOLD, TBORIG, TBBEST, TBBESTSCL)

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

END SUBROUTINE MOFIT
