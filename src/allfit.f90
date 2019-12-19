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

SUBROUTINE ALLFIT

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE NEBLISTARRAY
  USE UNIVARRAY
  USE PPOTARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, Z, NSNAP, ATMOL, TOTAT, COUNT, II, COUNT2, ACC
  INTEGER :: PICK, PARAMPICK, OLOOP, BACC
  INTEGER, ALLOCATABLE :: ATINMOL(:)
  REAL(LATTEPREC), ALLOCATABLE :: COORDS(:,:), FORCES(:,:), FORCEIN(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: BONDFORCES(:,:), FPAIR(:,:)
  REAL(LATTEPREC) :: Y, ENERGY, MYERR, MYERRINIT, MINERR, PREVERR
  REAL(LATTEPREC) :: RN, TMPERR
  REAL(LATTEPREC) :: BONDMINERR, BONDPREVERR
  REAL(LATTEPREC), ALLOCATABLE ::  PPBEST(:,:), PPOLD(:,:), PPKEEP(:,:)
  REAL(LATTEPREC), ALLOCATABLE ::  BINTBEST(:,:), BINTOLD(:,:)
  CHARACTER(LEN=2) :: X
  CHARACTER(LEN=2), ALLOCATABLE :: SPEC(:)
  IF (EXISTERROR) RETURN

  ALLOCATE(PPBEST(5,PP2FIT), PPOLD(5,PP2FIT), PPKEEP(5,PP2FIT))
  ALLOCATE(BINTBEST(2,BINT2FIT), BINTOLD(2,BINT2FIT))

  CALL INITRNG

  ! Read the xyz files and gradients

  OPEN(UNIT=60, STATUS="OLD", FILE="forces")

  NSNAP = PPNMOL*PPNGEOM

  TOTAT = 0
  DO I = 1, NSNAP
     !     print*, I
     READ(60,*) ATMOL
     TOTAT = TOTAT+ATMOL
     READ(60,*) ENERGY
     DO J  = 1, ATMOL
        READ(60,*) X, Y,Y,Y,Y,Y,Y
     ENDDO
  ENDDO

  REWIND(60)

  ALLOCATE(COORDS(3,TOTAT), FORCES(3,TOTAT), SPEC(TOTAT), ATINMOL(NSNAP))
  ALLOCATE(BONDFORCES(3,TOTAT), FPAIR(3,TOTAT), FORCEIN(3,TOTAT))

  COUNT = 0

  DO I = 1, NSNAP
     READ(60,*) ATINMOL(I)
     READ(60,*) ENERGY

     DO J = 1, ATINMOL(I)

        COUNT = COUNT + 1

        READ(60,*) SPEC(COUNT), COORDS(1,COUNT), COORDS(2,COUNT), &
             COORDS(3,COUNT), FORCEIN(1,COUNT), FORCEIN(2,COUNT), &
             FORCEIN(3,COUNT)

     ENDDO
  ENDDO

10 FORMAT(A1,6F12.6)

  FORCEIN = -FORCEIN*27.211385/0.591772

  ! First get the forces from the bond part and electrostatics

  BACC = 0

  DO OLOOP = 1, BINFITSTEP

     IF (OLOOP .GT. 1) THEN

        DO I = 1, BINT2FIT
           DO J = 1, 2
              BINTOLD(J,I) = BOND(J,I)
           ENDDO
        ENDDO

        CALL RANDOM_NUMBER(RN)

        PICK = INT(RN*REAL(BINT2FIT)) + 1

        CALL RANDOM_NUMBER(RN)

        PARAMPICK = INT(RN*TWO) + 1

        CALL RANDOM_NUMBER(RN)

        BOND(PARAMPICK,PICK) = BOND(PARAMPICK,PICK) * &
             (ONE + PPSIGMA*(TWO*RN - ONE))

        CALL UNIVTAILCOEF(BOND(:,PICK))

     ENDIF

     COUNT = 0
     COUNT2 = 0

     DO II = 1, NSNAP

        DEALLOCATE(CR, ATELE, F, FPP, FTOT)
        DEALLOCATE(DELTAQ, MYCHARGE)
        DEALLOCATE(ELEMPOINTER)
        IF (ELECTRO .EQ. 0) DEALLOCATE(LCNSHIFT)
        IF (KON .EQ. 1) DEALLOCATE(KF)
        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (SPINON .EQ. 0) THEN
              DEALLOCATE(FPUL, FSCOUL)
           ELSE
              DEALLOCATE(FPUL, FSCOUL, FSSPIN)
           ENDIF
        ENDIF

        ! Put the coordinates into CR

        NATS = ATINMOL(II)

        ALLOCATE(CR(3,NATS), ATELE(NATS), F(3,NATS), FPP(3,NATS), FTOT(3,NATS))
        ALLOCATE(DELTAQ(NATS), MYCHARGE(NATS))
        ALLOCATE(ELEMPOINTER(NATS))

        IF (ELECTRO .EQ. 0) THEN
           ALLOCATE(LCNSHIFT(NATS))
           LCNSHIFT = ZERO
        ENDIF

        IF (KON .EQ. 1) ALLOCATE(KF(3,NATS))

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (SPINON .EQ. 0) THEN
              ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS))
           ELSE
              ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS), FSSPIN(3,NATS))
           ENDIF
        ENDIF

        DO J = 1, NATS
           COUNT = COUNT + 1
           ATELE(J) = SPEC(COUNT)
           CR(1,J) = COORDS(1,COUNT)
           CR(2,J) = COORDS(2,COUNT)
           CR(3,J) = COORDS(3,COUNT)
        ENDDO

        ! Set up pointer to the data in TBparam/electrons.dat

        DO I = 1, NATS
           DO J = 1, NOELEM
              IF (ATELE(I) .EQ. ELE(J)) ELEMPOINTER(I) = J
           ENDDO
        ENDDO

        ! For use when getting the partial charges

        DEALLOCATE(QLIST)

        ! Allocate the Hamiltonian matrix

        IF (KON .EQ. 0) THEN

           ! Real space         
           DEALLOCATE(H, HDIAG)

        ELSE ! k-space  

           DEALLOCATE(HK, HKDIAG)

        ENDIF

        IF (BASISTYPE .EQ. "NONORTHO") DEALLOCATE(H0)

        IF (SPINON .EQ. 0) THEN

           ! No spins: allocate 1 double-occupied bond order matrix 

           IF (KON .EQ. 0) THEN

              DEALLOCATE(BO)

           ELSE

              DEALLOCATE(KBO)

           ENDIF

        ELSEIF (SPINON .EQ. 1) THEN

           DEALLOCATE(HUP, HDOWN)
           DEALLOCATE(RHOUP, RHODOWN)
           DEALLOCATE(H2VECT)

           IF (BASISTYPE .EQ. "NONORTHO") DEALLOCATE(SPINLIST)

           DEALLOCATE(DELTASPIN, OLDDELTASPIN)

        ENDIF

        CALL GETHDIM

        DEALLOCATE(MATINDLIST)

        IF (SPINON .EQ. 1) DEALLOCATE(SPININDLIST)

        CALL GETMATINDLIST

        IF (SPINON .EQ. 0) THEN
           DEALLOCATE(BOZERO)
        ELSE
           DEALLOCATE(RHOUPZERO, RHODOWNZERO)
        ENDIF

        CALL RHOZERO

        CALL GETBNDFIL


        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (II .EQ. 1) THEN
              CALL ALLOCATENONO
           ELSE
              CALL DEALLOCATENONO
              CALL ALLOCATENONO
           ENDIF
        ENDIF

        IF (II .EQ. 1) THEN

           CALL ALLOCATECOULOMB
           CALL INITCOULOMB

        ELSE 

           CALL DEALLOCATECOULOMB
           CALL ALLOCATECOULOMB
           CALL INITCOULOMB

        ENDIF

        IF (II .EQ. 1) THEN

           CALL ALLOCATENEBARRAYS           
           CALL NEBLISTS(0)

        ELSE 

           CALL DEALLOCATENEBARRAYS

           CALL ALLOCATENEBARRAYS

           DEALLOCATE(NEBTB, NEBPP)
           IF (ELECTRO .EQ. 1) DEALLOCATE(NEBCOUL)

           CALL NEBLISTS(0)

        ENDIF

        ! Now we have the arrays set up we can get the energy and forces

        ! Build the charge independent H matrix

        IF (KON .EQ. 0) THEN

           IF (SPONLY .EQ. 0) THEN
              CALL BLDNEWHS_SP
           ELSE
              CALL BLDNEWHS
           ENDIF

        ELSE

           CALL KBLDNEWH

        ENDIF

        IF (SPINON .EQ. 1) THEN
           CALL GETDELTASPIN
           CALL BLDSPINH
        ENDIF

        CALL ALLOCATEDIAG

        IF (ELECTRO .EQ. 0) CALL QNEUTRAL(0,1) ! Local charge neutrality

        IF (ELECTRO .EQ. 1) CALL QCONSISTENCY(0,1) ! Self-consistent charges     

        IF (KON .EQ. 0) THEN

           IF (SPONLY .EQ. 0) THEN
              CALL GRADHSP
           ELSE
              CALL GRADH
           ENDIF

        ELSE

           CALL KGRADH

        ENDIF

        FTOT = TWO * F

        IF (ELECTRO .EQ. 1) FTOT = FTOT + FCOUL

        IF (BASISTYPE .EQ. "NONORTHO") THEN

           CALL PULAY

           CALL FCOULNONO

           FTOT = FTOT - TWO*FPUL + FSCOUL

           IF (SPINON .EQ. 1) THEN
              CALL FSPINNONO
              FTOT = FTOT + FSSPIN
           ENDIF

        ENDIF

        DO I = 1, NATS
           COUNT2 = COUNT2 + 1
           BONDFORCES(1,COUNT2) = FTOT(1,I)
           BONDFORCES(2,COUNT2) = FTOT(2,I)
           BONDFORCES(3,COUNT2) = FTOT(3,I)
        ENDDO

        CALL DEALLOCATEDIAG                                                 

     ENDDO

     IF (BASISTYPE .EQ. "NONORTHO") CALL DEALLOCATENONO

     IF (ELECTRO .EQ. 1) CALL DEALLOCATECOULOMB

     CALL DEALLOCATENEBARRAYS

     DEALLOCATE(NEBTB, NEBPP)
     IF (ELECTRO .EQ. 1) DEALLOCATE(NEBCOUL)


     FORCES = FORCEIN - BONDFORCES


     ! Now for the optimization

     ! Initial error

     COUNT = 0
     COUNT2 = 0

     DO II = 1, NSNAP

        DEALLOCATE(CR, ATELE, FPP)

        NATS = ATINMOL(II)

        ALLOCATE(CR(3,NATS), ATELE(NATS), FPP(3,NATS))

        DO I = 1, NATS
           COUNT = COUNT+1
           CR(1,I) = COORDS(1,COUNT)
           CR(2,I) = COORDS(2,COUNT)
           CR(3,I) = COORDS(3,COUNT)
           ATELE(I) = SPEC(COUNT)
        ENDDO

        CALL PAIRPOTNONEB

        DO I = 1, NATS
           COUNT2 = COUNT2 + 1
           FPAIR(1,COUNT2) = FPP(1,I)
           FPAIR(2,COUNT2) = FPP(2,I)
           FPAIR(3,COUNT2) = FPP(3,I)
        ENDDO

     ENDDO

     MYERRINIT = ZERO
     DO I = 1, NSNAP

        TMPERR = ZERO
        DO J = 1, ATINMOL(I)

           TMPERR = TMPERR + &
                (FORCES(1,J) - FPAIR(1,J))*(FORCES(1,J) - FPAIR(1,J)) + &
                (FORCES(2,J) - FPAIR(2,J))*(FORCES(2,J) - FPAIR(2,J)) + &
                (FORCES(3,J) - FPAIR(3,J))*(FORCES(3,J) - FPAIR(3,J))

        ENDDO

        MYERRINIT = MYERRINIT + TMPERR/REAL(ATINMOL(I))

     ENDDO

     MYERRINIT = MYERRINIT/REAL(NSNAP)
     !  PRINT*, MYERRINIT

     !  ALLOCATE(PPORIG(5,PP2FIT), PPBEST(5,PP2FIT), PPOLD(5,PP2FIT))

     DO I = 1, PP2FIT
        DO J = 1, 5
           PPBEST(J,I) = POTCOEF(J,I)
        ENDDO
     ENDDO

     !     PPBEST = PPORIG

     ACC = 0
     DO Z = 1, PPNFITSTEP

        ! Change the pair potential and the cut-offs too

        DO I = 1, PP2FIT
           DO J = 1, 5
              PPOLD(J,I) = POTCOEF(J,I)
           ENDDO
        ENDDO

        CALL RANDOM_NUMBER(RN)

        PICK = INT(RN*REAL(PP2FIT)) + 1

        CALL RANDOM_NUMBER(RN)

        PARAMPICK = INT(RN*FIVE) + 1

        CALL RANDOM_NUMBER(RN)

        POTCOEF(PARAMPICK,PICK) = POTCOEF(PARAMPICK,PICK) * &
             (ONE + PPSIGMA*(TWO*RN - ONE))

        CALL VDWTAILCOEF

        ! Here we get the contribution to the forces from the
        ! pair potential

        COUNT = 0
        COUNT2 = 0

        IF (Z .EQ. 1) PREVERR = MYERRINIT
        IF (Z .EQ. 1) MINERR = MYERRINIT

        DO II = 1, NSNAP

           DEALLOCATE(CR, ATELE, FPP)

           NATS = ATINMOL(II)

           ALLOCATE(CR(3,NATS), ATELE(NATS), FPP(3,NATS))

           DO I = 1, NATS 
              COUNT = COUNT+1
              CR(1,I) = COORDS(1,COUNT)
              CR(2,I) = COORDS(2,COUNT)
              CR(3,I) = COORDS(3,COUNT)
              ATELE(I) = SPEC(COUNT)
           ENDDO

           CALL PAIRPOTNONEB

           DO I = 1, NATS
              COUNT2 = COUNT2 + 1
              FPAIR(1,COUNT2) = FPP(1,I)
              FPAIR(2,COUNT2) = FPP(2,I)
              FPAIR(3,COUNT2) = FPP(3,I)
           ENDDO

        ENDDO

        MYERR = ZERO
        DO I = 1, NSNAP

           TMPERR = ZERO
           DO J = 1, ATINMOL(I)

              TMPERR = TMPERR + &
                   (FORCES(1,J) - FPAIR(1,J))*(FORCES(1,J) - FPAIR(1,J)) + &
                   (FORCES(2,J) - FPAIR(2,J))*(FORCES(2,J) - FPAIR(2,J)) + &
                   (FORCES(3,J) - FPAIR(3,J))*(FORCES(3,J) - FPAIR(3,J))

           ENDDO

           MYERR = MYERR + TMPERR/REAL(ATINMOL(I))

        ENDDO

        MYERR = MYERR/REAL(NSNAP)

        !     PRINT*, MYERR

        IF (MYERR .LT. PREVERR) THEN

           ACC = ACC + 1
           PREVERR = MYERR

           !           PRINT*, ACC, MYERR

           IF (MYERR .LT. MINERR) THEN

              MINERR = MYERR

              DO I = 1, PP2FIT
                 DO J = 1, 5
                    PPBEST(J,I) = POTCOEF(J,I)
                 ENDDO
              ENDDO

           ENDIF

        ELSE

           CALL RANDOM_NUMBER(RN)

           IF (EXP(-(MYERR - PREVERR)*PPBETA) .GT. RN) THEN

              ACC = ACC+1
              PREVERR = MYERR

              !              PRINT*, ACC, MYERR

           ELSE

              ! Put the original coefficients back

              DO I = 1, PP2FIT
                 DO J = 1, 5
                    POTCOEF(J,I) = PPOLD(J,I)
                 ENDDO
              ENDDO

           ENDIF

        ENDIF

     ENDDO

     IF (OLOOP .EQ. 1) THEN

        BONDPREVERR = MINERR
        BONDMINERR = MINERR

        DO I = 1, BINT2FIT
           DO J = 1, 2
              BINTBEST(J,I) = BOND(J,I)
           ENDDO
        ENDDO

     ELSE

        IF (MINERR .LT. BONDPREVERR) THEN

           BACC = BACC + 1
           BONDPREVERR = MINERR

           PRINT*, BACC, MINERR

           IF (MINERR .LT. BONDMINERR) THEN

              BONDMINERR = MINERR

              DO I = 1, BINT2FIT
                 DO J = 1, 2
                    BINTBEST(J,I) = BOND(J,I)
                 ENDDO
              ENDDO

              DO I = 1, PP2FIT
                 DO J = 1, 5
                    PPKEEP(J,I) = POTCOEF(J,I)
                 ENDDO
              ENDDO

           ENDIF

        ELSE

           CALL RANDOM_NUMBER(RN)

           IF (EXP(-(MINERR - BONDPREVERR)*PPBETA) .GT. RN) THEN

              BACC = BACC+1
              BONDPREVERR = MINERR

              PRINT*, BACC, MINERR

           ELSE

              ! Put the original coefficients back

              DO I = 1, BINT2FIT
                 DO J = 1, 2
                    BOND(J,I) = BINTOLD(J,I)
                 ENDDO
              ENDDO

           ENDIF

        ENDIF

     ENDIF

  ENDDO

  DO I = 1, BINT2FIT
     WRITE(6,20) ELE1(I), ELE2(I), BINTBEST(1,I), BINTBEST(2,I)
  ENDDO


  DO I = 1, PP2FIT
     WRITE(6,20) PPELE1(I), PPELE2(I), PPKEEP(1,I), PPKEEP(2,I), &
          PPKEEP(3,I), PPKEEP(4,I), PPKEEP(5,I), POTCOEF(6,I), POTCOEF(7,I), &
          POTCOEF(8,I), POTCOEF(9,I), POTCOEF(10,I)

  ENDDO

  CALL ALLOCATENEBARRAYS

  CALL NEBLISTS(0)

  DEALLOCATE(PPKEEP, PPOLD, PPBEST)
  DEALLOCATE(BINTBEST, BINTOLD)

20 FORMAT(A2, 1X, A2, 1X, 10F16.8)

  RETURN

END SUBROUTINE ALLFIT
