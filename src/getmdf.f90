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

SUBROUTINE GETMDF(SWITCH, CURRITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE XBOARRAY
  USE MYPRECISION
  USE COULOMBARRAY
  USE DMARRAY, ONLY : HubForce, DELTAQDM
  USE TIMER_MOD
  USE MIXER_MOD

  IMPLICIT NONE

  INTEGER :: SWITCH, CURRITER, I, MLSI0, MLSI
  REAL(LATTEPREC) :: ZEROSCFMOD
  IF (EXISTERROR) RETURN

  MLSI0 = TIME_MLS()

  if (SWITCH .NE. 0) THEN

   IF(DOKERNEL .EQV. .TRUE.) THEN
   !IF(DOKERNEL .EQV. .TRUE. .AND. (.NOT.DFTBU))THEN 
      NRANK = MIN(NORECS,NATS)
      !!CALL KERNELPROPAGATION(CURRITER,NRANK)
      !IF ((CURRITER.GE.9).AND.(CURRITER.LT.12)) THEN
      !  WRITE(*,*) '# FULL K CURRITER = ', CURRITER
      !  CALL FULLKERNELPROPAGATION(CURRITER)
      !ELSEIF (MOD(CURRITER,200).EQ.1) THEN

      IF (KERNELSCHEME==1) THEN
        IF (CURRITER.LE.2) THEN
          WRITE(*,*) '# FULL K CURRITER = ', CURRITER
          CALL FULLKERNELPROPAGATION(CURRITER)
        ELSE
          WRITE(*,*) '# RANK_N K CURRITER = ', CURRITER
          CALL PRECONDKERNELPROPAGATION(CURRITER,NRANK)
        ENDIF
      ELSE  ! default
        ! new adaptive preconditioner kernel, Feb 29, 2020
        CALL ADAPTPRECONDKERNEL(CURRITER,NRANK)
      ENDIF
   ENDIF 

  ENDIF
 
  !
  ! The atoms have moved, so we have to build a new H (and overlap matrix)

  MLSI = TIME_MLS()
  
  IF (KON .EQ. 0) THEN
     IF(VERBOSE >= 1)WRITE(*,*)"KON = 0 ..."
     !IF (SPONLY .EQ. 0) THEN
     !   CALL BLDNEWHS_SP
     !ELSE
        CALL BLDNEWHS
     !ENDIF

  ELSE
     IF(VERBOSE >= 1)WRITE(*,*)"KON = 1 ..."
     CALL KBLDNEWH

  ENDIF
  FLUSH(6)

  !WRITE(*,*) "Time get H",  TIME_MLS() - MLSI

  MLSI = TIME_MLS()

  ! Broken?

  !  IF (SWITCH .EQ. 0 .AND. RESTART .EQ. 1) CALL IFRESTART

  IF (SPINON .EQ. 1 .AND. CURRITER .EQ. 1) THEN
     CALL GETDELTASPIN
     CALL BLDSPINH
  ENDIF

  IF (CONTROL .EQ. 5 .AND. CURRITER .EQ. 1) THEN

     ! We do this to initialize SP2 Fermi

     !     CONTROL = 2
     !     CALL QCONSISTENCY(0,1)
     !     CONTROL = 5
     CALL GERSHGORIN
     CALL SP2FERMIINIT
  ENDIF

  !  IF (ELECTRO .EQ. 1) THEN  ! Self-consistent charge transfer on

  !
  ! If we're running with XBO on, we start to propagate things
  ! after the first iteration
  !

  IF (XBOON .GT. 0 .AND. CURRITER .GT. 1 .AND. SWITCH .NE. 0) THEN

     IF(VERBOSE >= 1)WRITE(*,*)"XBOON .GT. 0 .AND. CURRITER .GT. 1 .AND. SWITCH .NE. 0 ..."
     IF (XBOON .EQ. 1) THEN ! We will add other types of XBO later

        ! Propagating partial charges or diagonal elements of H
        IF (DFTBU) THEN
          IF(VERBOSE >= 1)WRITE(*,*)"Doing XBODM ..."
          !CALL DMKERNELPROPAGATION(CURRITER)  !!! using rank-1 update
          CALL dP2MD(CURRITER)  !!! rank-m
          CALL XBODM(CURRITER) ! Propagate DM's
          DELTAQ = DELTAQDM    ! UPDATE DELTAQ FROM NEW DELTAQDM
        ELSE
          IF(VERBOSE >= 1)WRITE(*,*)"Doing XBO ..."
          CALL XBO(CURRITER) ! Propagate q's
        ENDIF
        !
        ! If we are also propagating the chemical potential
        !

        IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 .OR. CONTROL .EQ. 5) THEN

           CALL PROPCHEMPOT(CURRITER) ! Propagate mu

        ENDIF

        IF (SPINON .EQ. 1) CALL PROPSPINS(CURRITER) ! Propagate m's

        FLUSH(6)

     ENDIF

     FLUSH(6)
  ENDIF

  WRITE(*,*) "Time for PROPCHEMPOT XBO GETDELTASPIN",  TIME_MLS() - MLSI
  MLSI = TIME_MLS()
  !
  ! If SWITCH = 0, then we don't have a set of partials charges
  ! yet and we'll have to get them from the charge-independent
  ! H-matrix
  !
  ! Whether we're running full self consistency at each MD time step
  ! or only a user-specified number of iterations is determined by the
  ! value of FULLQCONV
  !

  IF (ELECTRO .EQ. 0) THEN
     IF(VERBOSE >= 1)WRITE(*,*)"Doing QNEUTRAL ..."
     CALL QNEUTRAL(SWITCH, CURRITER) ! Local charge neutrality
  ELSE
     IF(VERBOSE >= 1)WRITE(*,*)"Doing QCONSISTENCY ..."
     CALL QCONSISTENCY(SWITCH, CURRITER) ! Self consistent charge transfer
  ENDIF

  WRITE(*,*) "Time for QNEUTRAL QCONSISTENCY ",  TIME_MLS() - MLSI
  ! Run to self-consistency QITER = 0 -> only H(P) + D calculated 

  !
  ! Setting up our XBO arrays after the first iteration
  !

  ! We initialize once we have our first set of self-consistent q's, mu, and m's

  IF (XBOON .GT. 0 .AND. CURRITER .EQ. 1) THEN

     IF (XBOON .EQ. 1) THEN ! Other cases to come

        IF (DFTBU) THEN
          IF(VERBOSE >= 1)WRITE(*,*)"Doing XBODM ..."
          CALL XBODM(1)
        ELSE
          IF(VERBOSE >= 1)WRITE(*,*)"Doing XBO ..."
          CALL XBO(1)
        ENDIF

        IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 &
             .OR. CONTROL .EQ. 5) CALL PROPCHEMPOT(1)

        IF (SPINON .EQ. 1) CALL PROPSPINS(1)


     ENDIF

  ENDIF


  MLSI = TIME_MLS()
  ! Setting up the XBO arrays

  FLUSH(6)
  !
  ! Get the forces from the covalent part
  !

  ! When we're done with qconsistency we have the non-orthogonal density
  ! matrix so we don't have to de-orthogonalize again here

  IF(VERBOSE >= 1)WRITE(*,*)"Getting forces ..."

  IF (DFTBU) THEN
    call HUBBARDFORCE
  ENDIF

  CALL GETFORCE

  WRITE(*,*) "Time for GETFORCE  ",  TIME_MLS() - MLSI

  MLSI = TIME_MLS()

  IF (ELECTRO .EQ. 1 .AND. QITER .EQ. 0) THEN

     OLDDELTAQS = DELTAQ ! save the propagated charges
     CALL GETDELTAQ ! Get updated set of partial charges

     ECOUL = ZERO

     FTOT = FTOT - FCOUL ! We're going to correct the electrostatic force

     DO I = 1, NATS

        ZEROSCFMOD = (TWO*DELTAQ(I) - OLDDELTAQS(I))/OLDDELTAQS(I)
        FCOUL(1,I) = FCOUL(1,I)*ZEROSCFMOD
        FCOUL(2,I) = FCOUL(2,I)*ZEROSCFMOD
        FCOUL(3,I) = FCOUL(3,I)*ZEROSCFMOD
        IF (DFTBU) THEN
          ECOUL = ECOUL + (TWO*DELTAQ(I) - OLDDELTAQS(I)) * &    ! Works only with EBand0 energy from TRRHOH0 = Tr[D*H0]
               (HUBBARDU(ELEMPOINTER(I))*OLDDELTAQS(I) + COULOMBV(I))
        ELSE
          ECOUL = ECOUL + OLDDELTAQS(I)* &    ! Orig, ecp Works with regular Eband and TRRHOH = tr[D*H]
              (HUBBARDU(ELEMPOINTER(I))*OLDDELTAQS(I) + COULOMBV(I))
        ENDIF

     ENDDO

     ECOUL = ECOUL/TWO

     FTOT = FTOT + FCOUL

  ENDIF

  WRITE(*,*)"Time for GETMDF", TIME_MLS()-MLSI0

  FLUSH(6)

  RETURN

END SUBROUTINE GETMDF
