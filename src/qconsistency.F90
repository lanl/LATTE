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

SUBROUTINE QCONSISTENCY(SWITCH, MDITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE COULOMBARRAY
  USE TIMER_MOD
  USE MYPRECISION
  USE MIXER_MOD
  USE DMARRAY ! CHANGE ANDERS
  USE NONOARRAY ! CHANGE ANDERS TEMP FOR CHECK ONLY

  IMPLICIT NONE

  INTEGER :: I, SWITCH, MDITER, ITER, II
  INTEGER :: ALLOKQ, ALLOKM, ALLOK, MDSOFT
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX, ITERACC
  REAL(4) :: TIMEACC
  REAL(LATTEPREC) :: MAXDQ
  REAL(LATTEPREC), ALLOCATABLE :: QDIFF(:), SPINDIFF(:)

  MDSOFT = 1   ! CHANGE ANDERS

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

  TIMEACC = 0.0
  ITERACC = 0
  write(*,*) 'QCONSISTENCY  MDITER = ', MDITER
  write(*,*) 'QCONSISTENCY  MDSOFT = ', MDSOFT

  IF (FULLQCONV .EQ. 1 .OR. MDITER .LE. MDSOFT) THEN

     IF (SWITCH .EQ. 0) THEN

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN

#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#elif defined(PROGRESSOFF)

              CALL ORTHOMYH
#endif

           ELSEIF (KON .EQ. 1) THEN
              CALL KORTHOMYH
           ENDIF
        ENDIF

        write(*,*) ' QCONS ORTHOH = ', ORTHOH(1,1:HDIM)

!        DOrth_old = BO  ! ANDERS CHANGE


        ! Compute the density matrix

        TX = START_TIMER(DMBUILD_TIMER)

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        DOrth_old = BO  ! ANDERS CHANGE
        write(*,*) ' QCONS DORTH_Old = ', BO(1,1:HDIM)

        TX = STOP_TIMER(DMBUILD_TIMER)

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
#elif defined(PROGRESSOFF)
              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF

        IF(VERBOSE >= 1)WRITE(*,*)"In GETDELTAQ ..."
        CALL GETDELTAQ

        IF (SPINON .EQ. 1) CALL GETDELTASPIN

     ENDIF

     !
     ! Now we're going to run our iterations for self-consistency
     !

     ALLOK = 1
     ITER = 0

     ALLOCATE(QDIFF(NATS))
     IF (SPINON .EQ. 1) ALLOCATE(SPINDIFF(DELTADIM))

     DO WHILE (ALLOK .GT. 0)

        ITER = ITER + 1
        !IF(VERBOSE >= 1)WRITE(*,*)"SCF_ITER =", ITER   !! ANDERS
        FLUSH(6)


        IF (ELECMETH .EQ. 0) THEN

           ! First do the real space part of the electrostatics
           ! This subroutine is based on Sanville's work

           CALL COULOMBRSPACE

           ! And now the long range bit (this is also a modified version
           ! of Ed's code).

           CALL COULOMBEWALD

        ELSE

           ! Doing the electrostatics all in real space
           ! help tremendously when getting the virial

           CALL GASPCOULOMB

        ENDIF

        ! Now let's modify the diagonal elements of our H matrix according
        ! to the electrostatic potential experienced by each atom

        CALL ADDQDEP

        write(*,*) 'innan rad 1 H = ', H(1,1:HDIM)  ! ANDERS CHANGE
        write(*,*) 'innan rad 2 H = ', H(2,1:HDIM)  ! ANDERS CHANGE
        CALL ADDDFTBU_INIT ! ANDERS CHANGE
        write(*,*) ' rad 1 H = ', H(1,1:HDIM)  ! ANDERS CHANGE
        write(*,*) ' rad 2 H = ', H(2,1:HDIM)  ! ANDERS CHANGE

        ! Got to add the electrostatic potential to
        ! the Slater-Koster H before adding H_2 to form
        ! H_up and H_down

        !
        ! Calculate the spin-dependent H matrix again
        !

        IF (SPINON .EQ. 1) CALL BLDSPINH


        ! We've made changes to the H matrix so we have to re-orthogonalize

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#elif defined(PROGRESSOFF)
              CALL ORTHOMYH
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KORTHOMYH
           ENDIF
        ENDIF
        write(*,*) ' HOrth = ', ORTHOH(1,1:HDIM)  ! ANDERS CHANGE
        write(*,*) ' SMAT = ', SMAT(1,1:HDIM)  ! ANDERS CHANGE
        write(*,*) ' XMAT = ', XMAT(1,1:HDIM)  ! ANDERS CHANGE

        !
        ! New Hamiltonian: get the bond order
        !

        ! Compute the density matrix

        TX = START_TIMER(DMBUILD_TIMER)

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        DOrth = BO  ! ANDERS CHANGE
        write(*,*) ' QConsist BO orth = ', BO(1,1:HDIM)

        TX = STOP_TIMER(DMBUILD_TIMER)

        ! We have a density matrix computed in from the orthogonalized
        ! H matrix - we have to revert back

        !
        ! Save our old charges/spins so we can mix them later
        !

!        CALL ENTROPY

        IF (BASISTYPE.EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
                 CALL DEORTHOMYRHOPRG
              ELSE
                 CALL DEORTHOMYRHO
              ENDIF
#elif defined(PROGRESSOFF)

              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF

        OLDDELTAQS = DELTAQ

        !
        ! Get a new set of charges for our system
        !

        CALL GETDELTAQ

        !
        ! Let's check for convergence
        !

        ALLOKQ = 0

        QDIFF = ABS(DELTAQ - OLDDELTAQS)

        MAXDQ = MAXVAL(QDIFF)
        IF (MAXDQ .GT. ELEC_QTOL) ALLOKQ = 1

        ! Mix new and old partial charges

        IF (MDITER .LE. MDSOFT) THEN
#ifdef PROGRESSON
           IF(MX%MIXERON)THEN
              CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
           ELSE
              DELTAQ = QMIX*DELTAQ + (ONE - QMIX)*OLDDELTAQS
           ENDIF
#elif defined(PROGRESSOFF)
          ! NRANK = MIN(NORECS,NATS)         !ANDERS: STEALING UNUSED NORECS INPUT FOR NRANK NEEDS FIX!
          ! NRANK = 0 -> linear mixing, NRANK = NATS gives exact Newton
          ! CALL KERNELMIXER(ITER,NRANK)     !ANDERS: Alternative mixing scheme ANDERS
          ! DELTAQ = QMIX*DELTAQ + (ONE - QMIX)*OLDDELTAQS
                     !!! ANDERS CHECK

           write(*,*) ' QCONSISTENCY MIX A DOrth_old = ', DOrth_old(1,1:HDIM)
           write(*,*) ' QCONSISTENCY MIX A Delta_DOrth = ', DOrth(1,1:HDIM) - DOrth_old(1,1:HDIM)
           BO = DOrth_old + QMIX*(DOrth - DOrth_old)  ! ANDERS CHANGE
           DOrth_old = BO                             ! ANDERS CHANGE
           CALL DEORTHOMYRHO
           CALL GETDELTAQ

           !!! ANDERS CHECK
!           write(*,*) ' TBMD HUBBARD ENERGY in QCONS'
!           call HUBBARDFORCE
!           write(*,*) ' TBMD HUBBARD ENERGY = ', EHub

#endif

        ELSE

#ifdef PROGRESSON
           IF(MX%MIXERON)THEN
              CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
           ELSE
              DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
           ENDIF
#elif defined(PROGRESSOFF)
           DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
           BO = DOrth_old + QMIX*(DOrth - DOrth_old)  ! ANDERS CHANGE
           DOrth_old = BO                             ! ANDERS CHANGE
           CALL DEORTHOMYRHO
!          CALL GETDELTAQ
#endif

        ENDIF

        DOrth = DOrth_old  ! ANDERS CHANGE

        IF(VERBOSE >= 1) WRITE(*,*)"SCF_Iter:", ITER,", SCF error (MAXDQ) =",MAXDQ

        ALLOKM = 0

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           CALL GETDELTASPIN

           SPINDIFF = ABS(DELTASPIN - OLDDELTASPIN)

           IF (MAXVAL(SPINDIFF) .GT. SPINTOL)  ALLOKM = 1

           ! Mix new and old spin densities

           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN

        ENDIF

        ALLOK = ALLOKQ + ALLOKM

!        MAXSCF = 3  ! CHANGE ANDERS

        IF (ITER .EQ. MAXSCF) THEN

           WRITE(6,*) "# WARNING - the SCF procedure has not converged"
           WRITE(6,*) "# to the tolerances defined in TBparam/control.in"
           WRITE(6,*) "# Continuing anyway, but be very careful... "

           ALLOK = 0

           IF (STOPATMAXSCF) CALL ERRORS("qconsistency","The SCF procedure has not converged")

        ENDIF


     ENDDO

     NUMSCF = ITER

     DEALLOCATE(QDIFF)
     IF (SPINON .EQ. 1) DEALLOCATE(SPINDIFF)

     !     WRITE(6,99) "# HDIM, TIME PER RHO BUILD = ", HDIM, TIMEACC/REAL(ITERACC)
     !99   FORMAT(A29, I5, 1X, G18.6)

     FLUSH(6)

  ELSEIF (FULLQCONV .EQ. 0 .AND. MDON .EQ. 1 .AND. MDITER .GT. MDSOFT) THEN

     ! Now we're doing MD

     DO II = 1, QITER

        IF (ELECMETH .EQ. 0) THEN

           CALL COULOMBRSPACE

           CALL COULOMBEWALD

        ELSE

           CALL GASPCOULOMB

        ENDIF

        CALL ADDQDEP
        CALL ADDDFTBU  ! ANDERS CHANGE

        !
        ! Building the spin up and spin down H's after we've
        ! added the electrostatic potential to the Slater-Koster one,
        ! as it should be.
        !

        IF (SPINON .EQ. 1) CALL BLDSPINH

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN

#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
                 CALL ORTHOMYHPRG
              ELSE
                 CALL ORTHOMYH
              ENDIF
#elif defined(PROGRESSOFF)
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

        OLDDELTAQS = DELTAQ
        DOrth = BO  ! ANDERS CHANGE

        !
        ! Get a new set of charges for our system
        !

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
              IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
                 CALL DEORTHOMYRHOPRG
              ELSE
                 CALL DEORTHOMYRHO
              ENDIF
#elif defined(PROGRESSOFF)
              CALL DEORTHOMYRHO
#endif
           ELSEIF (KON .EQ. 1) THEN
              CALL KDEORTHOMYRHO
           ENDIF
        ENDIF

        CALL GETDELTAQ

        !
        ! Mix to get new charges
        !

        ! Linear mixing BOMD with QITER >= 1
        DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS

        BO = DOrth_old + QMIX*(DOrth - DOrth_old)  ! ANDERS CHANGE
        DOrth_old = BO                             ! ANDERS CHANGE
        CALL DEORTHOMYRHO
!       CALL GETDELTAQ

        !        PRINT*, DELTAQ(1)

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           CALL GETDELTASPIN
           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN

        ENDIF

     ENDDO
     ! Here we start if QITER = 0 in MD

     ! Calculate the bond order one more time since we need the forces for
     ! that charge distribution

     IF (ELECMETH .EQ. 0) THEN

        CALL COULOMBRSPACE

        CALL COULOMBEWALD

     ELSE

        CALL GASPCOULOMB

     ENDIF

     CALL ADDQDEP
     CALL ADDDFTBU  ! ANDERS CHANGE

     ! This is the right order

     IF (SPINON .EQ. 1) CALL BLDSPINH

     IF (BASISTYPE .EQ. "NONORTHO") THEN
        IF (KON .EQ. 0) THEN

#ifdef PROGRESSON
           IF (LATTEINEXISTS) THEN  !orthogonalize from progress lib if latte.in exists
              CALL ORTHOMYHPRG
           ELSE
              CALL ORTHOMYH
           ENDIF
#elif defined(PROGRESSOFF)
           CALL ORTHOMYH
#endif


        ELSEIF (KON .EQ. 1) THEN
           CALL KORTHOMYH
        ENDIF
     ENDIF


     !
     ! New Hamiltonian: get the bond order/density matrices
     !

     ! Compute the density matrix

     IF (KON .EQ. 0) THEN
        CALL GETRHO(MDITER)
     ELSE
        CALL KGETRHO
     ENDIF

     DOrth = BO

     write(*,*) 'QCONSISTENCYCHECK- DOING_MD? --'  ! ANDERS

     IF (BASISTYPE .EQ. "NONORTHO") THEN
        IF (KON .EQ. 0) THEN
#ifdef PROGRESSON
           IF (LATTEINEXISTS) THEN  !deorthogonalize from progress lib if latte.in exists
              CALL DEORTHOMYRHOPRG
           ELSE
              CALL DEORTHOMYRHO
           ENDIF
#elif defined(PROGRESSOFF)
           CALL DEORTHOMYRHO
#endif
        ELSEIF (KON .EQ. 1) THEN
           CALL KDEORTHOMYRHO
        ENDIF
     ENDIF

  ENDIF

  RETURN

END SUBROUTINE QCONSISTENCY
