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

SUBROUTINE QNEUTRAL(SWITCH, MDITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, SWITCH, MDITER, ITER, II
  INTEGER :: ALLOKQ, ALLOKM, ALLOK
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
  REAL(4) :: TIMEACC
  REAL(LATTEPREC), ALLOCATABLE :: SPINDIFF(:)
  IF (EXISTERROR) RETURN

  !
  ! If FULLQCONV = 1, then we're going to iterate until all charges are within
  ! QTOL.
  !
  ! If FULLQCONV = 0, then we're going to run only a user specified number
  ! of iterations (= QITER)
  !
  ! If SWITCH = 0, then we don't have any partial charges defined yet so
  ! we'll have to get these from our charge-independent H matrix first
  !

  ENTE = ZERO

  IF (FULLQCONV .EQ. 1 .OR. MDITER .LE. 10) THEN

     IF (SWITCH .EQ. 0) THEN

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#else
              CALL ORTHOMYH
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KORTHOMYH
           ENDIF
        ENDIF

        ! Compute the density matrix

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        ! If we used diagonalization, we can compute the response function
        ! for updating the charges

        !        IF (CONTROL .EQ. 1) CALL GETRESPF


        !
        ! Now we have our bond-order/density matrices,
        ! we can get the charges and spins
        !

        ! The partial charges and spin densities are computed from the
        ! density matrices computed from the orthogonalized H matrices

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
                 CALL DEORTHOMYRHOPRG
              ELSE
                 CALL DEORTHOMYRHO
              ENDIF
#else
              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF

        CALL GETDELTAQ

        IF (SPINON .EQ. 1) CALL GETDELTASPIN

     ENDIF

     !
     ! Now we're going to run our iterations for self-consistency
     !

     ALLOK = 1
     ITER = 0

     IF (SPINON .EQ. 1) ALLOCATE(SPINDIFF(DELTADIM))

     DO WHILE (ALLOK .GT. 0)

        ITER = ITER + 1

        ! Adjust on-site energies such that partial charges -> 0

        CALL SHIFTH(QMIX)

        !
        ! New Hamiltonian: get the bond order
        !

        ! Compute the density matrix

        IF (SPINON .EQ. 1) CALL BLDSPINH

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#else
              CALL ORTHOMYH
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KORTHOMYH
           ENDIF
        ENDIF


        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        !        IF (CONTROL .EQ. 1) CALL GETRESPF

        IF (BASISTYPE.EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
                 CALL DEORTHOMYRHOPRG
              ELSE
                 CALL DEORTHOMYRHO
              ENDIF
#else
              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF


        CALL GETDELTAQ

        IF (ABS(MAXVAL(DELTAQ)) .LT. ELEC_QTOL) ALLOK = 0
        !PRINT*, ABS(MAXVAL(DELTAQ)), ELEC_QTOL

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           CALL GETDELTASPIN

           SPINDIFF = ABS(DELTASPIN - OLDDELTASPIN)

           IF (MAXVAL(SPINDIFF) .LT. SPINTOL)  ALLOK = 0

           ! Mix new and old spin densities

           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN

        ENDIF


        IF (ITER .EQ. MAXSCF) THEN

           WRITE(6,*) "WARNING - the SCF procedure has not converged"
           WRITE(6,*) "to the tolerances defined in TBparam/control.in"
           WRITE(6,*) "Continuing anyway, but be very careful... "

           ALLOK = 0

           IF (STOPATMAXSCF) CALL ERRORS("qneutral","The SCF procedure has not converged")

        ENDIF


     ENDDO

     NUMSCF = ITER

     IF (SPINON .EQ. 1) DEALLOCATE(SPINDIFF)

  ELSEIF (FULLQCONV .EQ. 0 .AND. MDON .EQ. 1 .AND. MDITER .GT. 10) THEN

     ! Now we're doing MD

     DO II = 1, QITER

        CALL SHIFTH(MDMIX)

        IF (SPINON .EQ. 1) CALL BLDSPINH

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#else
              CALL ORTHOMYH
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KORTHOMYH
           ENDIF
        ENDIF



        !
        ! New Hamiltonian: get the bond order
        !

        ! Compute the density matrix

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
                 CALL DEORTHOMYRHOPRG
              ELSE
                 CALL DEORTHOMYRHO
              ENDIF
#else
              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF

        !        IF (CONTROL .EQ. 1) CALL GETRESPF

        CALL GETDELTAQ

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           CALL GETDELTASPIN
           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN

        ENDIF

     ENDDO

  ENDIF

  RETURN

END SUBROUTINE QNEUTRAL
