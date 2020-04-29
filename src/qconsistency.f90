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
  USE DMARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, SWITCH, MDITER, ITER, II, DIIS_DIM, INFO
  INTEGER :: ALLOKQ, ALLOKM, ALLOK, NEW_MIXER, MDSOFT
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX, ITERACC, DIIS_M
  INTEGER, ALLOCATABLE :: IPIV(:)
  REAL(4) :: TIMEACC
  REAL(LATTEPREC) :: MAXDQ, MAG_DIFF, DELTA_DIFF, SCF_ERR, MLSI
  REAL(LATTEPREC), ALLOCATABLE :: QDIFF(:), SPINDIFF(:)
  REAL(LATTEPREC), ALLOCATABLE :: QDIFF_DIIS(:,:), BMAT_DIIS(:,:), DIIS_RHS(:)
  REAL(LATTEPREC), ALLOCATABLE :: QHIST(:,:)

  MDSOFT = 10
  IF (DFTBU) THEN
    MDSOFT = 1 
    SCF_ERR = 1.D0
  ENDIF

  IF (EXISTERROR) RETURN

  ALLOCATE(QDIFF_DIIS(NATS, MAXSCF), QHIST(NATS, MAXSCF))

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

  NEW_MIXER = 0

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

        IF (DFTBU) DOrth_old = BO

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

     DIIS_M = 0

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

        IF (DFTBU) CALL ADDDFTBU(.true.) 

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

        IF (DFTBU) DOrth = BO 

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

        IF (ITER .GT. 1) THEN
           DELTA_DIFF = MAG_DIFF
        ELSE
           DELTA_DIFF = 1.0D6
        ENDIF

        MAG_DIFF = ZERO
        DO I = 1, NATS
           MAG_DIFF = MAG_DIFF + QDIFF(I)*QDIFF(I)
        ENDDO

        MAG_DIFF = SQRT(MAG_DIFF/REAL(NATS))
        
        DELTA_DIFF = DELTA_DIFF - MAG_DIFF

        ALLOKQ = 0
        DO I = 1, NATS
           IF (ABS(QDIFF(I)) .GT. ELEC_QTOL) ALLOKQ = ALLOKQ + 1
        ENDDO

        IF ((DELTA_DIFF .GT. 1.0D-3) .AND. (NEW_MIXER .EQ. 0)) THEN

           ! Run linear mixing for a bit until we start to converge

        !IF (MDITER .LE. MDSOFT) THEN
          IF (.NOT.DFTBU) THEN
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
           !CALL DMKERNELMIXER(ITER,MAXDQ)    ! Rank-1 updated DM kernel mixer
           CALL dP2MIXER(ITER,SCF_ERR,MAXDQ)  ! Rank-m updated DM kernel mixer
           DOrth_old = BO                   
           CALL DEORTHOMYRHO
           OLDDELTAQS = DELTAQ
           CALL GETDELTAQ
           SCF_ERR = norm2(DELTAQ - OLDDELTAQS)/sqrt(ONE*NATS)
          ENDIF

        ELSE

           ! Flip this flag so we don't go back into linear mixing

           NEW_MIXER = 1

           IF (DFTBU) THEN
             BO = DOrth_old + QMIX*(DOrth - DOrth_old)
             DOrth_old = BO
             CALL DEORTHOMYRHO
             CALL GETDELTAQ      ! INCLUDED_GETDELTAQ Probably not to comment out
           ELSE
   
           ! Then switch to DIIS

           DIIS_M = DIIS_M + 1


           ! Linear mixing for the first iteration in DIIS

           IF (DIIS_M .EQ. 1) THEN

              ! The history of the errors

              QDIFF_DIIS(:,DIIS_M) = DELTAQ - OLDDELTAQS

              ! The history of the deltaq's. We use DIIS to combine these into a better guess

              QHIST(:,DIIS_M) = DELTAQ

              DELTAQ = QMIX*DELTAQ + (ONE - QMIX)*OLDDELTAQS
           
           ELSE

              QDIFF_DIIS(:,DIIS_M) = DELTAQ - OLDDELTAQS
              QHIST(:,DIIS_M) = DELTAQ

              ALLOCATE(BMAT_DIIS(DIIS_M+1, DIIS_M+1), DIIS_RHS(DIIS_M+1), IPIV(DIIS_M+1))
              
              ! Setting up 'B'

              BMAT_DIIS = MINUSONE
              BMAT_DIIS(DIIS_M+1, DIIS_M+1) = ZERO
              
              DO I = 1, DIIS_M
                 DO J = 1, DIIS_M
                    
                    DO K = 1, NATS
                       
                       BMAT_DIIS(J,I) = BMAT_DIIS(J,I) + QDIFF_DIIS(K,J)*QDIFF_DIIS(K,I)
                       
                    ENDDO
                    
                 ENDDO
              ENDDO
              
              DIIS_DIM = DIIS_M + 1
              
              DIIS_RHS = ZERO
              DIIS_RHS(DIIS_M+1) = MINUSONE

              ! Solve A X = B using LAPACK call
              
              CALL DGESV(DIIS_DIM, 1, BMAT_DIIS, DIIS_DIM, IPIV, DIIS_RHS, DIIS_DIM, INFO)
              
!              PRINT*, INFO, DIIS_M, DIIS_RHS, DELTAQ(1)
              
              DELTAQ = ZERO
              DO I = 1, DIIS_M
                 DELTAQ(:) = DELTAQ(:) + QHIST(:,I)*DIIS_RHS(I)
              ENDDO
              
              DEALLOCATE(BMAT_DIIS, DIIS_RHS, IPIV)
              
           ENDIF
           ENDIF ! DFTBU
        ENDIF

        IF (DFTBU)  DOrth = DOrth_old

        IF(VERBOSE >= 1)WRITE(*,*)"SCF error (MAXDQ) =",MAXDQ

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

        IF (ITER .EQ. MAXSCF) THEN

           WRITE(6,*) "# WARNING - the SCF procedure has not converged"
           WRITE(6,*) "# to the tolerances defined in TBparam/control.in"
           WRITE(6,*) "# Continuing anyway, but be very careful... "

           ALLOK = 0

           IF (STOPATMAXSCF) CALL ERRORS("qconsistency","The SCF procedure has not converged")

!           DO I = 1, NATS
!              PRINT*, ATELE(I), CR(1,I), CR(2,I), CR(3,I)
!           ENDDO
!           PRINT*, CHARGE
           
        ENDIF

        IF (MDON .EQ. 1 .AND. ITER .LE. QITER) ALLOK = 1


     ENDDO

     NUMSCF = ITER

     DEALLOCATE(QDIFF)
     IF (SPINON .EQ. 1) DEALLOCATE(SPINDIFF)

     !     WRITE(6,99) "# HDIM, TIME PER RHO BUILD = ", HDIM, TIMEACC/REAL(ITERACC)
     !99   FORMAT(A29, I5, 1X, G18.6)

     IF (MDON .EQ. 1 .AND. MDADAPT .EQ. 1) THEN

        IF (EGAP .GT. MINGAP) THEN

           MDADAPT_COUNT = MDADAPT_COUNT + 1
!           PRINT*, MDADAPT_COUNT, EGAP
           IF (MDADAPT_COUNT .EQ. 100) THEN

              WRITE(6,'("# Returning to linear mixing")')
              FULLQCONV = 0
              QITER = 1

           ENDIF


        ENDIF

     ENDIF

     FLUSH(6)

  ELSEIF (FULLQCONV .EQ. 0 .AND. MDON .EQ. 1 .AND. MDITER .GT. MDSOFT) THEN

     IF (MDON .EQ. 1 .AND. MDADAPT .EQ. 1) THEN

        ! When the gap closes start running to full convergence                     

        IF (EGAP .LT. MINGAP) THEN
        
           WRITE(6,'("# Gap closed! Running to full SCF convergence")')
           
           FULLQCONV = 1
           QITER = 4
!           MDADAPT = 0
           MDADAPT_COUNT = 0

        ENDIF


     ENDIF

     ! Now we're doing MD

     DO II = 1, QITER

        IF (ELECMETH .EQ. 0) THEN

           CALL COULOMBRSPACE

           CALL COULOMBEWALD

        ELSE

           CALL GASPCOULOMB

        ENDIF

        CALL ADDQDEP
        IF (DFTBU) CALL ADDDFTBU(.false.) 

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
        IF (DFTBU) DOrth = BO  

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

        DELTAQ = MDMIX*DELTAQ + (ONE - MDMIX)*OLDDELTAQS

        ! Here we do DM mixing instead of charge mixing
        IF (DFTBU) THEN
          BO = DOrth_old + QMIX*(DOrth - DOrth_old) 
          DOrth_old = BO                           
          CALL DEORTHOMYRHO
          CALL GETDELTAQ         ! INCLUDED_GETDELTAQ
        ENDIF

        !        PRINT*, DELTAQ(1)

        IF (SPINON .EQ. 1) THEN

           OLDDELTASPIN = DELTASPIN
           CALL GETDELTASPIN
           DELTASPIN = SPINMIX*DELTASPIN + (ONE - SPINMIX)*OLDDELTASPIN

        ENDIF

     ENDDO

     ! Calculate the bond order one more time since we need the forces for
     ! that charge distribution

     IF (ELECMETH .EQ. 0) THEN

        CALL COULOMBRSPACE

        CALL COULOMBEWALD

     ELSE

        CALL GASPCOULOMB

     ENDIF

     CALL ADDQDEP
     IF (DFTBU) CALL ADDDFTBU(.false.) 

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

     IF (DFTBU) THEN
       DOrth_old = DOrth
       DOrth = BO
     ENDIF

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

  DEALLOCATE(QDIFF_DIIS, QHIST)

  RETURN

END SUBROUTINE QCONSISTENCY
