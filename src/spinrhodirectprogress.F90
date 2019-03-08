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

SUBROUTINE SPINRHOEVECSPRG

#ifdef PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE SPINARRAY
  USE MDARRAY
  USE MYPRECISION
  USE SPARSEARRAY

  USE BML

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: dp = LATTEPREC
  INTEGER :: I, J, K
  INTEGER :: ITER, BREAKLOOP
  INTEGER :: UPNE, DOWNNE
  REAL(LATTEPREC) :: OCCERROR , OCCUP, OCCDOWN
  REAL(LATTEPREC) :: EXPARG, FDIRAC
  REAL(LATTEPREC) :: SHIFTCP
  REAL(LATTEPREC) :: FDIRACARG, DFDIRAC
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: HOMO, LUMO
  REAL(LATTEPREC), ALLOCATABLE :: FDIRACUP(:), FDIRACDOWN(:)
  REAL(LATTEPREC) :: S, OCCLOGOCC_ELECTRONS, OCCLOGOCC_HOLES
  TYPE(BML_MATRIX_T) :: RHO_BML, AUX_BML, AUX1_BML, EVECS_BML, OCCUPATION_BML

  IF (EXISTERROR) RETURN

  !  NUMLIMIT = EXP(-EXPTOL)

  !RHOUP = ZERO
  !RHODOWN = ZERO

  IF (VERBOSE >= 1) WRITE(*,*)"In SPINRHOEVECSPRG ..."

  ALLOCATE(FDIRACUP(HDIM))
  ALLOCATE(FDIRACDOWN(HDIM))

  ITER = 0

  OCCERROR = 100000000.0

  !
  ! For a finite electronic temperature
  !

  IF (KBT .GT. 0.000001) THEN

     DO WHILE ( ABS(OCCERROR) .GT. BREAKTOL .AND. ITER .LT. 100 )

        ITER = ITER + 1

        OCCUP = ZERO
        OCCDOWN = ZERO
        DFDIRAC = ZERO

        DO I = 1, HDIM

           FDIRACARG = (UPEVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           EXPARG = EXP(FDIRACARG)
           FDIRAC = ONE/(ONE + EXPARG)
           OCCUP = OCCUP + FDIRAC
           DFDIRAC = DFDIRAC + EXPARG*FDIRAC*FDIRAC

           FDIRACARG = (DOWNEVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)


           EXPARG = EXP(FDIRACARG)
           FDIRAC = ONE/(ONE + EXPARG)
           OCCDOWN = OCCDOWN + FDIRAC
           DFDIRAC = DFDIRAC + EXPARG*FDIRAC*FDIRAC

        ENDDO

        DFDIRAC = DFDIRAC/KBT

        OCCERROR = TOTNE - OCCUP - OCCDOWN

        IF (ABS(DFDIRAC) .LT. NUMLIMIT) DFDIRAC = SIGN(NUMLIMIT, DFDIRAC)

        SHIFTCP = OCCERROR/DFDIRAC

        IF (ABS(SHIFTCP) .GT. MAXSHIFT) SHIFTCP = SIGN(MAXSHIFT, SHIFTCP)

        CHEMPOT = CHEMPOT + SHIFTCP

     ENDDO

     IF (ITER .GT. 100) THEN
        CALL ERRORS("spinrhodirect","Newton-Raphson scheme to find the &
             & chemical potential has not converged")
     ENDIF


     ! Now we have the chemical potential we can build the density matrices

     ! Entropy and HOMO-LUMO: only computed when we need them during MD

     S = ZERO

     IF (MDON .EQ. 0 .OR. &
          (MDON .EQ. 1 .AND. MOD(ENTROPYITER, WRTFREQ) .EQ. 0 )) THEN

        ! We will compute the HOMO-LUMO gap and the entropy in here too

        ! Starting guesses

        HOMO = MIN(UPEVALS(1), DOWNEVALS(1))
        LUMO = MAX(UPEVALS(HDIM), DOWNEVALS(HDIM))

        DO I = 1, HDIM

           ! Entropy

           FDIRACARG = (UPEVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           FDIRAC = ONE/(ONE + EXP(FDIRACARG))

           OCCLOGOCC_ELECTRONS = FDIRAC * LOG(FDIRAC)
           OCCLOGOCC_HOLES = (ONE - FDIRAC) * LOG(ONE - FDIRAC)

           S = S + OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES

           FDIRACARG = (DOWNEVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           FDIRAC = ONE/(ONE + EXP(FDIRACARG))

           OCCLOGOCC_ELECTRONS = FDIRAC * LOG(FDIRAC)
           OCCLOGOCC_HOLES = (ONE - FDIRAC) * LOG(ONE - FDIRAC)

           S = S + OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES

           ! Homo-Lumo

           IF (UPEVALS(I) .LT. CHEMPOT .AND. &
                UPEVALS(I) .GT. HOMO) HOMO = UPEVALS(I)
           IF (DOWNEVALS(I) .LT. CHEMPOT .AND. &
                DOWNEVALS(I) .GT. HOMO) HOMO = DOWNEVALS(I)

           ! LUMO = smallest eigenvalue > chemical potential

           IF (UPEVALS(I) .GT. CHEMPOT .AND. &
                UPEVALS(I) .LT. LUMO) LUMO = UPEVALS(I)
           IF (DOWNEVALS(I) .GT. CHEMPOT .AND. &
                DOWNEVALS(I) .LT. LUMO) LUMO = DOWNEVALS(I)


        ENDDO

        EGAP = LUMO - HOMO

     ENDIF

     ENTE = -KBT*S

     ! Building diagonal

     DO I = 1, HDIM

        FDIRACARG = (UPEVALS(I) - CHEMPOT)/KBT

        FDIRACARG = MAX(FDIRACARG, -EXPTOL)
        FDIRACARG = MIN(FDIRACARG, EXPTOL)

        FDIRACUP(I) = ONE/(ONE + EXP(FDIRACARG))

        FDIRACARG = (DOWNEVALS(I) - CHEMPOT)/KBT

        FDIRACARG = MAX(FDIRACARG, -EXPTOL)
        FDIRACARG = MIN(FDIRACARG, EXPTOL)

        FDIRACDOWN(I) = ONE/(ONE + EXP(FDIRACARG))

     ENDDO
     

  ELSE  ! For zero electronic temperature

     ! Occupy the lowest eigenvalues

     ! See whether we occupy spin up or down first

     IF (UPEVALS(1) .LT. DOWNEVALS(1)) THEN
        UPNE = 1
        DOWNNE = 0
        CHEMPOT = UPEVALS(1)
     ELSE
        UPNE = 0
        DOWNNE = 1
        CHEMPOT = DOWNEVALS(1)
     ENDIF

     ! Go through the eigenvalues and occupy if
     ! 1) its energy >= the previous occupied eigenvalue
     ! 2) its energy <= the energy of the unoccupied eigenvalue from the other
     ! spin channel

     ! Stop when the total number of electrons reaches the target
     ! Works with O and O2.

     BREAKLOOP = 0

     DO WHILE ( BREAKLOOP .EQ. 0 )

        IF (UPEVALS(UPNE + 1) .LE. DOWNEVALS(DOWNNE + 1)) THEN
           UPNE = UPNE + 1
           CHEMPOT = UPEVALS(UPNE)
           FDIRACUP(UPNE) = 1.0_dp 
        ENDIF

        IF ( UPNE + DOWNNE .EQ. INT(TOTNE) ) BREAKLOOP = 1

        IF (BREAKLOOP .EQ. 0 .AND. &
             DOWNEVALS(DOWNNE + 1) .LE. UPEVALS(UPNE + 1)) THEN
           DOWNNE = DOWNNE + 1
           CHEMPOT = DOWNEVALS(DOWNNE)
           FDIRACDOWN(DOWNNE) = 1.0_dp 
        ENDIF

        IF ( UPNE + DOWNNE .EQ. INT(TOTNE) ) BREAKLOOP = 1

!        print*, upne, downne, totne

     ENDDO

     ! Here HOMO = CHEMPOT, so the energy of the
     ! LUMO = min(upevals(upne + 1), downevals(downne + 1)

     LUMO = MIN(UPEVALS(UPNE + 1), DOWNEVALS(DOWNNE + 1))

     EGAP = LUMO - CHEMPOT

  ENDIF

  !Rho up
  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, RHO_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, OCCUPATION_BML)
  CALL BML_SET_DIAGONAL(OCCUPATION_BML,FDIRACUP)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, AUX_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, EVECS_BML)
  CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
            UPEVECS, EVECS_BML, ZERO, HDIM)

  CALL BML_MULTIPLY(EVECS_BML,OCCUPATION_BML,AUX_BML,1.0_dp,0.0_dp,NUMTHRESH)
  CALL BML_DEALLOCATE(OCCUPATION_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, AUX1_BML)
  CALL BML_TRANSPOSE(EVECS_BML, AUX1_BML)
  CALL BML_DEALLOCATE(EVECS_BML) 

  CALL BML_MULTIPLY(AUX_BML,AUX1_BML,RHO_BML,1.0_dp,0.0_dp,NUMTHRESH)
  CALL BML_DEALLOCATE(AUX_BML) 
  CALL BML_DEALLOCATE(AUX1_BML) 

  CALL BML_EXPORT_TO_DENSE(RHO_BML, RHOUP)
  CALL BML_DEALLOCATE(RHO_BML)

  !Rho down
  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, RHO_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, OCCUPATION_BML)
  CALL BML_SET_DIAGONAL(OCCUPATION_BML,FDIRACDOWN)


  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, AUX_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, EVECS_BML)
  CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
            DOWNEVECS, EVECS_BML, ZERO, HDIM)

  CALL BML_MULTIPLY(EVECS_BML,OCCUPATION_BML,AUX_BML,1.0_dp,0.0_dp,NUMTHRESH)
  CALL BML_DEALLOCATE(OCCUPATION_BML)

  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, AUX1_BML)
  CALL BML_TRANSPOSE(EVECS_BML, AUX1_BML)
  CALL BML_DEALLOCATE(EVECS_BML)

  CALL BML_MULTIPLY(AUX_BML,AUX1_BML,RHO_BML,1.0_dp,0.0_dp,NUMTHRESH)
  CALL BML_DEALLOCATE(AUX_BML)
  CALL BML_DEALLOCATE(AUX1_BML)

  CALL BML_EXPORT_TO_DENSE(RHO_BML, RHODOWN)
  CALL BML_DEALLOCATE(RHO_BML)


  DEALLOCATE(FDIRACUP)
  DEALLOCATE(FDIRACDOWN)

  RETURN

#endif

END SUBROUTINE SPINRHOEVECSPRG
