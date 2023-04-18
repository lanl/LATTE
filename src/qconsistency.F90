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
  USE KSPACEARRAY, ONLY : KBO
  USE DMARRAY
  USE LATTEPARSER 
#ifdef  PROGRESSON
  USE SPARSEARRAY, ONLY : NUMTHRESH
  USE BML
  USE NONOARRAYPROGRESS
#endif

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = LATTEPREC
  INTEGER :: I, SWITCH, MDITER, ITER, II
  INTEGER :: ALLOKQ, ALLOKM, ALLOK
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX, ITERACC
  REAL(4) :: TIMEACC
  LOGICAL :: FIRSTCALL=.TRUE.
  REAL(LATTEPREC) :: MAXDQ, SCF_ERR, MLSI, MLSII, MAXDSP, MLSI0
  REAL(LATTEPREC), ALLOCATABLE :: QDIFF(:), SPINDIFF(:)

  MDSOFT = 10
  IF (DFTBU) THEN
    MDSOFT = 1 
    SCF_ERR = 1.D0
  ENDIF

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

        ! Compute the density matrix

        TX = START_TIMER(DMBUILD_TIMER)

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

#ifdef PROGRESSON
        IF (DFTBU .AND. KON==0) CALL BML_COPY_NEW(ORTHOBO_BML,DO_BML_OLD)
#elif defined(PROGRESSOFF)
        IF (DFTBU .AND. KON==0) DOrth_old = BO
#endif

        !IF (DFTBU .AND. KON==1) DORK_OLD = KBO

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
        IF(VERBOSE >= 1)WRITE(*,*)"SCF ITER =", ITER
        FLUSH(6)


        IF (ELECMETH .EQ. 0) THEN

           !
           ! First do the real space part of the electrostatics
           ! This subroutine is based on Sanville's work
           !

           CALL COULOMBRSPACE

           !
           ! And now the long range bit (this is also a modified version
           ! of Ed's code).
           !

           CALL COULOMBEWALD

        ELSE

           !
           ! Doing the electrostatics all in real space
           ! help tremendously when getting the virial
           !

           CALL GASPCOULOMB

        ENDIF

        !
        ! Now let's modify the diagonal elements of our H matrix according
        ! to the electrostatic potential experienced by each atom
        !

        CALL ADDQDEP

#ifdef PROGRESSON
        IF (DFTBU) CALL ADDDFTBUPRG(.true.) 
#elif defined(PROGRESSOFF)
        IF (DFTBU) CALL ADDDFTBU(.true.) 
#endif

        ! Got to add the electrostatic potential to
        ! the Slater-Koster H before adding H_2 to form
        ! H_up and H_down

        !
        ! Calculate the spin-dependent H matrix again
        !

        IF (SPINON .EQ. 1) CALL BLDSPINH


        ! We've made changes to the H matrix so we have to re-orthogonalize
        MLSI = TIME_MLS()

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

        IF(VERBOSE >= 1) WRITE(*,*) "Time for GETMDF-QCONSISTENCY-ORTHOMYH ",  TIME_MLS() - MLSI
        !
        ! New Hamiltonian: get the bond order
        !

        ! Compute the density matrix

        MLSI = TIME_MLS()
        TX = START_TIMER(DMBUILD_TIMER)

        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        IF (DFTBU .AND. KON==0) then 
#ifdef  PROGRESSON
           !CALL BML_COPY_NEW(ORTHOBO_BML,PNO_BML)
           CALL BML_COPY_NEW(ORTHOBO_BML,DO_BML)
#elif defined(PROGRESSOFF)
           DOrth = BO 
           !PNO   = BO
#endif
        ENDIF

        !IF (DFTBU .AND. KON==1) DORK = KBO

        IF(VERBOSE >=1) WRITE(*,*) "Time for GETMDF-QCONSISTENCY-BUILDRHO ",  TIME_MLS() - MLSI
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
!           IF (DFTBU .AND. KON==0) PNO = BO
!#ifdef  PROGRESSON
!           IF (DFTBU .AND. KON==0) CALL BML_COPY_NEW(BO_BML,PNO_BML)
!#endif
        ENDIF

        OLDDELTAQS = DELTAQ

        !
        ! Get a new set of charges for our system
        !

        !CALL GETDELTAQ
        IF (.NOT.DFTBU) CALL GETDELTAQ

        !
        ! Let's check for convergence
        !

        ALLOKQ = 0

        IF (.NOT.DFTBU) THEN
          QDIFF = ABS(DELTAQ - OLDDELTAQS)

          MAXDQ = MAXVAL(QDIFF)
        ENDIF

        !IF (MAXDQ .GT. ELEC_QTOL) ALLOKQ = 1

        ! Mix new and old partial charges

        IF (MDITER .LE. MDSOFT) THEN
          IF (.NOT.DFTBU .OR. KON==1) THEN
           ! mixing charge
           !IF (.NOT.DFTBU) THEN
#ifdef PROGRESSON
           IF(MX%MIXERON)THEN
              CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
           ELSE
              DELTAQ = QMIX*DELTAQ + (ONE - QMIX)*OLDDELTAQS
           ENDIF
#elif defined(PROGRESSOFF)
           DELTAQ = QMIX*DELTAQ + (ONE - QMIX)*OLDDELTAQS
#endif
          ELSE
           !IF (KON==0) THEN
             IF (DOKERNEL) THEN
               !CALL DMKERNELMIXER(ITER,MAXDQ)              ! Rank1-1 updated DM kernel mixer
               CALL dP2MIXER(ITER,SCF_ERR,MAXDQ)            ! Rank1-1 updated DM kernel mixer
             ELSE
#ifdef PROGRESSON
               IF(MX%MIXERON)THEN
                 ! use pulayDM
                 CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
               ELSE
                 CALL BML_ADD(DO_BML,DO_BML_OLD,1.0_DP,-1.0_DP,NUMTHRESH)
                 CALL BML_ADD(DO_BML,DO_BML_OLD,QMIX,1.0_DP,NUMTHRESH)
                 CALL BML_COPY_NEW(DO_BML, ORTHOBO_BML)
               ENDIF
#elif defined(PROGRESSOFF)
               DOrth = DOrth_old + QMIX*(DOrth - DOrth_old)   ! Simple linear mixing
               BO = DOrth
#endif
             ENDIF

#ifdef PROGRESSON
             CALL BML_COPY_NEW(ORTHOBO_BML, DO_BML_OLD)
             CALL DEORTHOMYRHOPRG
#elif defined(PROGRESSOFF)
             DOrth_old = BO
             CALL DEORTHOMYRHO
#endif

             OLDDELTAQS = DELTAQ
             CALL GETDELTAQ
             SCF_ERR = norm2(DELTAQ - OLDDELTAQS)/sqrt(ONE*NATS) ! ANDERS CHANGE SCF_ERR TO COMPARE TO Dev Verison
             MAXDQ = MAXVAL(abs(DELTAQ-OLDDELTAQS))

           !ELSE
           !  KBO = DORK_old + QMIX*(DORK - DORK_old)   ! Simple linear mixing
           !  DORK_old = KBO                            !
           !  CALL KDEORTHOMYRHO
           !  OLDDELTAQS = DELTAQ
           !  CALL GETDELTAQ
           !ENDIF
          ENDIF
        ELSE

          IF (.NOT.DFTBU .OR. KON==1) THEN
#ifdef PROGRESSON
           IF(MX%MIXERON)THEN
              CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
           ELSE
              DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
           ENDIF
#elif defined(PROGRESSOFF)
           DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
#endif
          ELSE
             ! mixing DM
#ifdef PROGRESSON
             IF(MX%MIXERON)THEN
               ! use pulayDM
               CALL QMIXPRG(ITER)     !Alternative mixing scheme from PROGRESS
             ELSE
               CALL BML_ADD(DO_BML,DO_BML_OLD,1.0_DP,-1.0_DP,NUMTHRESH)
               CALL BML_ADD(DO_BML,DO_BML_OLD,QMIX,1.0_DP,NUMTHRESH)
               CALL BML_COPY_NEW(DO_BML, ORTHOBO_BML)
             ENDIF
#elif defined(PROGRESSOFF)
             DOrth = DOrth_old + QMIX*(DOrth - DOrth_old)   ! Simple linear mixing
             BO = DOrth
#endif

#ifdef PROGRESSON
             CALL BML_COPY_NEW(ORTHOBO_BML, DO_BML_OLD)
             CALL DEORTHOMYRHOPRG
#elif defined(PROGRESSOFF)
             DOrth_old = BO
             CALL DEORTHOMYRHO
#endif
             OLDDELTAQS = DELTAQ
             CALL GETDELTAQ      ! INCLUDED_GETDELTAQ Probably not to comment out ANDERS?
             MAXDQ = MAXVAL(abs(DELTAQ-OLDDELTAQS))
          ENDIF

        ENDIF

#ifdef PROGRESSON
        IF (DFTBU.AND.KON==0) CALL BML_COPY_NEW(DO_BML_OLD, DO_BML)
#elif defined(PROGRESSOFF)
        IF (DFTBU.AND.KON==0) DOrth = DOrth_old
#endif

        IF(VERBOSE >= 1)WRITE(*,*)"SCF error (MAXDQ) =",MAXDQ," SCF Tol =",ELEC_QTOL

        IF (MAXDQ .GT. ELEC_QTOL) ALLOKQ = 1

        ALLOKM = 0

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           OLDSUMSPIN   = SUMSPIN

           CALL GETDELTASPIN

           SPINDIFF = ABS(DELTASPIN - OLDDELTASPIN)
           MAXDSP = MAXVAL(SPINDIFF)

           IF (MAXDSP .GT. SPINTOL)  ALLOKM = 1

           ! Mix new and old spin densities

           IF (MAXDSP .GT. 1.0D-1 .OR. ITER < 3 .OR. (.NOT. DOKERNEL)) THEN
             DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN
             SUMSPIN   = SPINMIX*SUMSPIN   + (ONE - SPINMIX)*OLDSUMSPIN
             FIRSTCALL = .TRUE.
           ELSE
             CALL ADAPTKERNELMIXER(FIRSTCALL, NORECS)
             FIRSTCALL = .FALSE.
           ENDIF
           IF (DOKERNEL)  CALL REDUCE_DELTASPIN(NATS, DELTADIM, SUMSPIN, DELTAQ, 2)

           IF(VERBOSE >= 1)WRITE(*,*)"SPIN error =",MAXDSP," SPIN Tol =",SPINTOL

        ENDIF

        ALLOK = ALLOKQ + ALLOKM

        IF (ITER .EQ. MAXSCF) THEN

           WRITE(6,*) "# WARNING - the SCF procedure has not converged"
           WRITE(6,*) "# to the tolerances defined in TBparam/control.in"
           WRITE(6,*) "# Continuing anyway, but be very careful... "

           ALLOK = 0

           IF (STOPATMAXSCF) CALL ERRORS("qconsistency","The SCF procedure has not converged")

        ENDIF


     ENDDO

     IF (DFTBU .AND. KON==0) PNO = BO
#ifdef  PROGRESSON
     IF (DFTBU .AND. KON==0) CALL BML_COPY_NEW(BO_BML,PNO_BML)
#endif

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

           MLSII = TIME_MLS()
           CALL COULOMBRSPACE
           WRITE(*,*)"Time for COULOMBRSPACE",TIME_MLS()-MLSII

           MLSII = TIME_MLS()
           CALL COULOMBEWALD
           WRITE(*,*)"Time for COULOMBEWALD",TIME_MLS()-MLSII

        ELSE

           CALL GASPCOULOMB

        ENDIF

        CALL ADDQDEP
#ifdef PROGRESSON
        IF (DFTBU) CALL ADDDFTBUPRG(.false.) 
#elif defined(PROGRESSOFF)
        IF (DFTBU) CALL ADDDFTBU(.false.) 
#endif

        !
        ! Building the spin up and spin down H's after we've
        ! added the electrostatic potential to the Slater-Koster one,
        ! as it should be.
        !

        IF (SPINON .EQ. 1) CALL BLDSPINH

        MLSI = TIME_MLS()
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

        IF (VERBOSE >=1 )WRITE(*,*)"Time for GETMDF-QCONSISTENCY-ORTHOMYH", TIME_MLS() - MLSI

        !
        ! New Hamiltonian: get the bond order
        !

        ! Compute the density matrix


        MLSI = TIME_MLS()
        IF (KON .EQ. 0) THEN
           CALL GETRHO(MDITER)
        ELSE
           CALL KGETRHO
        ENDIF

        IF (VERBOSE>=1) WRITE(*,*)"Time for GETMDF-QCONSISTENCY-GETRHO", TIME_MLS() - MLSI
        OLDDELTAQS = DELTAQ
#ifdef PROGRESSON
        IF (DFTBU .AND. KON==0) CALL BML_COPY_NEW(ORTHOBO_BML,DO_BML)
#elif defined(PROGRESSOFF)
        IF (DFTBU .AND. KON==0) DOrth = BO
#endif
        !IF (DFTBU .AND. KON==1) DORK = KBO

        !
        ! Get a new set of charges for our system
        !

        MLSI = TIME_MLS()
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
        IF (VERBOSE >=1 ) WRITE(*,*)"Time for GETMDF-QCONSISTENCY-DEORTHOMYRHO", TIME_MLS() - MLSI

        CALL GETDELTAQ

        !
        ! Mix to get new charges
        !

#ifdef PROGRESSON
        IF(MX%MIXERON)THEN
          CALL QMIXPRG(II)
        ELSE
          DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
        ENDIF
#else
        DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS
#endif

        ! Here we do DM mixing instead of charge mixing
        IF (DFTBU .AND. KON==0) THEN
#ifdef PROGRESSON
          CALL BML_COPY(DO_BML_OLD, ORTHOBO_BML)
          CALL BML_ADD(ORTHOBO_BML,DO_BML_OLD,1.0_DP,-1.0_DP,NUMTHRESH)
          CALL BML_ADD(ORTHOBO_BML,DO_BML_OLD,QMIX,1.0_DP,NUMTHRESH)
          CALL BML_COPY_NEW(ORTHOBO_BML, DO_BML_OLD)
          CALL DEORTHOMYRHOPRG
#elif defined(PROGRESSOFF)
          CALL dP2MIXER(ITER,SCF_ERR,MAXDQ)            ! Rank1-1 updated DM kernel mixer
          !BO = DOrth_old + QMIX*(DOrth - DOrth_old) 
          DOrth_old = BO                           
          CALL DEORTHOMYRHO
#endif
          CALL GETDELTAQ         ! INCLUDED_GETDELTAQ
        ENDIF
        !IF (DFTBU .AND. KON==1) THEN
        !  KBO = DORK_old + QMIX*(DORK - DORK_old) 
        !  DORK_old = KBO                           
        !  CALL DEORTHOMYRHO
        !  CALL GETDELTAQ         ! INCLUDED_GETDELTAQ
        !ENDIF

        !        PRINT*, DELTAQ(1)

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           OLDSUMSPIN   = SUMSPIN

           CALL GETDELTASPIN
           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN
           SUMSPIN   = SPINMIX*SUMSPIN   + (ONE - SPINMIX)*OLDSUMSPIN

           IF (DOKERNEL)  CALL REDUCE_DELTASPIN(NATS,DELTADIM,SUMSPIN, DELTAQ,2)
        ENDIF

     ENDDO

     ! Calculate the bond order one more time since we need the forces for
     ! that charge distribution

     IF (ELECMETH .EQ. 0) THEN

        MLSII = TIME_MLS()
        CALL COULOMBRSPACE
        WRITE(*,*)"Time for COULOMBRSPACE",TIME_MLS()-MLSII

        MLSII = TIME_MLS()
        CALL COULOMBEWALD
        WRITE(*,*)"Time for COULOMBEWALD",TIME_MLS()-MLSII


     ELSE

        CALL GASPCOULOMB

     ENDIF

     CALL ADDQDEP
#ifdef PROGRESSON
     IF (DFTBU) CALL ADDDFTBUPRG(.false.) 
     !IF (DFTBU) CALL ADDDFTBU(.false.) 
#elif defined(PROGRESSOFF)
     IF (DFTBU) CALL ADDDFTBU(.false.) 
#endif

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

     IF (DFTBU .AND. KON==0) THEN
#ifdef PROGRESSON
       CALL BML_COPY_NEW(DO_BML, DO_BML_OLD)
       CALL BML_COPY_NEW(ORTHOBO_BML, DO_BML)
#elif defined(PROGRESSOFF)
       DOrth_old = DOrth
       DOrth = BO
#endif
     ENDIF

     !IF (DFTBU .AND. KON==1) THEN
     !  DORK_old = DORK
     !  DORK = KBO
     !ENDIF

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
