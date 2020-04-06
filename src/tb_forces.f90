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

SUBROUTINE TBFORCES

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE NONOARRAY
  USE NEBLISTARRAY
  USE SPINARRAY
  USE VIRIALARRAY
  USE MYPRECISION
  USE DMARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, L, M, N, KK, INDI, INDJ
  INTEGER :: LBRA, MBRA, LKET, MKET
  INTEGER :: PREVJ, NEWJ
  INTEGER :: PBCI, PBCJ, PBCK
  INTEGER :: BASISI(5), BASISJ(5), LBRAINC, LKETINC
  INTEGER :: SPININDI, SPININDJ
  REAL(LATTEPREC) :: ALPHA, BETA, PHI, RHO, RHODIFF, COSBETA
  REAL(LATTEPREC) :: RIJ(3), DC(3)
  REAL(LATTEPREC) :: MAGR, MAGR2, MAGRP, MAGRP2
  REAL(LATTEPREC) :: FTMP_BOND(3)
  REAL(LATTEPREC) :: FTMP_PULAY(3), FTMP_COUL(3), FTMP_SPIN(3)
  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2
  REAL(LATTEPREC) :: WSPINI, WSPINJ
  REAL(LATTEPREC) :: MYDFDA, MYDFDB, MYDFDR, RCUTTB
  REAL(LATTEPREC) :: SMYDFDA, SMYDFDB, SMYDFDR
  LOGICAL PATH
  IF (EXISTERROR) RETURN

  F = ZERO
  FPUL = ZERO
  VIRPUL = ZERO
  VIRBOND = ZERO

  ! We first have to make the matrix S^-1 H rho = X^2 H rho

  IF (SPINON .EQ. 0) THEN

     IF (KBT .GT. 0.0000001) THEN

#ifdef DOUBLEPREC

#ifdef GPUON

        ! XMAT * XMAT = S^-1                                                                                              

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, XMAT, XMAT, X2HRHO)

        ! S^-1 * H                                                                                                        

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, X2HRHO, H, NONOTMP)

        ! (S^-1 * H)*RHO                                                                                                  

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, NONOTMP, BO, X2HRHO)

#elif defined(GPUOFF)

        ! XMAT * XMAT = S^-1

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             XMAT, HDIM, XMAT, HDIM, ZERO, X2HRHO, HDIM)

        ! S^-1 * H

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             X2HRHO, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)

        ! (S^-1 * H)*RHO

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, BO, HDIM, ZERO, X2HRHO, HDIM)

#endif
        
#elif defined(SINGLEPREC)

        ! XMAT * XMAT = S^-1

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             XMAT, HDIM, XMAT, HDIM, ZERO, X2HRHO, HDIM)

        ! S^-1 * H

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             X2HRHO, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)

        ! (S^-1 * H)*RHO

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, BO, HDIM, ZERO, X2HRHO, HDIM)

#endif

     ELSE

        ! Te = 0 : Fp = 2Tr[rho H rho dS/dR]

        ! Be careful - we're working with bo = 2rho, so we need 
        ! the factor of 1/2...

#ifdef DOUBLEPREC

#ifdef GPUON

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, BO, H, NONOTMP)

        CALL mmlatte(HDIM, 0, 0, HALF, ZERO, NONOTMP, BO, X2HRHO)

#elif defined(GPUOFF)


        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             BO, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, HALF, &
             NONOTMP, HDIM, BO, HDIM, ZERO, X2HRHO, HDIM)

#endif

#elif defined(SINGLEPREC)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             BO, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, HALF, &
             NONOTMP, HDIM, BO, HDIM, ZERO, X2HRHO, HDIM)

#endif        

     ENDIF

  ELSE

     IF (KBT .GT. 0.0000001) THEN

#ifdef DOUBLEPREC

#ifdef GPUON

        ! XMAT * XMAT = S^-1                                                                                              

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, XMAT, XMAT, SPINTMP)

        !      Hup * rhoup                                                                                                

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, HUP, RHOUP, NONOTMP)

        ! (Hup * rhoup) + Hdown*rhodown                                                                                   

        CALL mmlatte(HDIM, 0, 0, ONE, ONE, HDOWN, RHODOWN, NONOTMP)

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, SPINTMP, NONOTMP, X2HRHO)

#elif defined(GPUOFF)


        ! XMAT * XMAT = S^-1

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             XMAT, HDIM, XMAT, HDIM, ZERO, SPINTMP, HDIM)

        !      Hup * rhoup

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             HUP, HDIM, RHOUP, HDIM, ZERO, NONOTMP, HDIM)

        ! (Hup * rhoup) + Hdown*rhodown 

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             HDOWN, HDIM, RHODOWN, HDIM, ONE, NONOTMP, HDIM)


        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             SPINTMP, HDIM, NONOTMP, HDIM, ZERO, X2HRHO, HDIM)

#endif

#elif defined(SINGLEPREC)

        ! XMAT * XMAT = S^-1

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             XMAT, HDIM, XMAT, HDIM, ZERO, SPINTMP, HDIM)

        !      Hup * rhoup

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             HUP, HDIM, RHOUP, HDIM, ZERO, NONOTMP, HDIM)

        ! (Hup * rhoup) + Hdown*rhodown 

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             HDOWN, HDIM, RHODOWN, HDIM, ONE, NONOTMP, HDIM)


        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             SPINTMP, HDIM, NONOTMP, HDIM, ZERO, X2HRHO, HDIM)

#endif

     ELSE

#ifdef DOUBLEPREC

        ! At zero K: tr(rho H rho)

#ifdef GPUON

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, RHOUP, HUP, NONOTMP)

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, NONOTMP, RHOUP, X2HRHO)

        CALL mmlatte(HDIM, 0, 0, ONE, ZERO, RHODOWN, HDOWN, NONOTMP)

        CALL mmlatte(HDIM, 0, 0, ONE, ONE, NONOTMP, RHODOWN, X2HRHO)

#elif defined(GPUOFF)

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             RHOUP, HDIM, HUP, HDIM, ZERO, NONOTMP, HDIM)

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, RHOUP, HDIM, ZERO, X2HRHO, HDIM)

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             RHODOWN, HDIM, HDOWN, HDIM, ZERO, NONOTMP, HDIM)

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, RHODOWN, HDIM, ONE, X2HRHO, HDIM)

#endif

#elif defined(SINGLEPREC)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             RHOUP, HDIM, HUP, HDIM, ZERO, NONOTMP, HDIM)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, RHOUP, HDIM, ZERO, X2HRHO, HDIM)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             RHODOWN, HDIM, HDOWN, HDIM, ZERO, NONOTMP, HDIM)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
             NONOTMP, HDIM, RHODOWN, HDIM, ONE, X2HRHO, HDIM)

#endif

     ENDIF

  ENDIF

  IF (DFTBU) X2HRHO = X2HRHO - 0.5D0 * AHub

!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &
!$OMP SHARED(CR, BOX, X2HRHO) &
!$OMP SHARED(BO, RHOUP, RHODOWN, SPINON, ELECTRO) & 
!$OMP SHARED(HCUT, SCUT, MATINDLIST, BASISTYPE) &
!$OMP SHARED(DELTASPIN, WSS, WPP, WDD, WFF, SPININDLIST, H2VECT) &
!$OMP SHARED(HUBBARDU, DELTAQ, COULOMBV, ORBITAL_LIST, CUTOFF_LIST) & 
!$OMP SHARED(LCNSHIFT) & 
!$OMP SHARED(FPUL, F) &
!$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &
!$OMP PRIVATE(RIJ, MAGR2, MAGR, MAGRP2, MAGRP, PATH, PHI, ALPHA, BETA, COSBETA) &
!$OMP PRIVATE(FTMP_PULAY, FTMP_COUL, FTMP_SPIN, FTMP_BOND) &
!$OMP PRIVATE(DC, LBRAINC, LBRA, MBRA, L, LKETINC, LKET, MKET, RHO, RHODIFF) &
!$OMP PRIVATE(MYDFDA, MYDFDB, MYDFDR, RCUTTB) &
!$OMP PRIVATE(SMYDFDA, SMYDFDB, SMYDFDR)&
!$OMP REDUCTION(+:VIRPUL, VIRBOND)


  DO I = 1, NATS

     ! Build list of orbitals on atom I

     BASISI(:) = ORBITAL_LIST(:,I)

     INDI = MATINDLIST(I)
!     IF (SPINON .EQ. 1) SPININDI = SPININDLIST(I)

     DO NEWJ = 1, TOTNEBTB(I)

        J = NEBTB(1, NEWJ, I)
        PBCI = NEBTB(2, NEWJ, I)
        PBCJ = NEBTB(3, NEWJ, I)
        PBCK = NEBTB(4, NEWJ, I)        

        RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1) - CR(1,I)

        RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2) - CR(2,I)

        RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3) - CR(3,I)

        MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

        ! Building the forces is expensive - use the cut-off

        RCUTTB = CUTOFF_LIST(J,I)

        IF (MAGR2 .LT. RCUTTB*RCUTTB) THEN

           MAGR = SQRT(MAGR2)


           BASISJ(:) = ORBITAL_LIST(:,J)

           INDJ = MATINDLIST(J)
!           IF (SPINON .EQ. 1) SPININDJ = SPININDLIST(J)


           MAGRP2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2)
           MAGRP = SQRT(MAGRP2)

           ! transform to system in which z-axis is aligned with RIJ

           PATH = .FALSE.
           IF (ABS(RIJ(1)) .GT. 1E-12) THEN
              IF (RIJ(1) .GT. ZERO .AND. RIJ(2) .GE. ZERO) THEN
                 PHI = ZERO
              ELSEIF (RIJ(1) .GT. ZERO .AND. RIJ(2) .LT. ZERO) THEN
                 PHI = TWO * PI
              ELSE
                 PHI = PI
              ENDIF
              ALPHA = ATAN(RIJ(2) / RIJ(1)) + PHI
           ELSEIF (ABS(RIJ(2)) .GT. 1E-12) THEN
              IF (RIJ(2) .GT. 1E-12) THEN
                 ALPHA = PI / TWO
              ELSE
                 ALPHA = THREE * PI / TWO
              ENDIF
           ELSE
              ! pathological case: pole in alpha at beta=0
              PATH = .TRUE.
           ENDIF

           COSBETA = RIJ(3)/MAGR
           BETA = ACOS(RIJ(3) / MAGR) 

           DC = RIJ/MAGR

           ! build forces using PRB 72 165107 eq. (12) - the sign of the
           ! dfda contribution seems to be wrong, but gives the right 
           ! answer(?)

           FTMP_BOND = ZERO
           FTMP_PULAY = ZERO
           FTMP_COUL = ZERO
           FTMP_SPIN = ZERO

           K = INDI

           LBRAINC = 1
           DO WHILE (BASISI(LBRAINC) .NE. -1)

              LBRA = BASISI(LBRAINC)
              LBRAINC = LBRAINC + 1

              DO MBRA = -LBRA, LBRA

                 K = K + 1
                 L = INDJ

                 LKETINC = 1
                 DO WHILE (BASISJ(LKETINC) .NE. -1)

                    LKET = BASISJ(LKETINC)
                    LKETINC = LKETINC + 1

                    DO MKET = -LKET, LKET

                       L = L + 1

                       SELECT CASE(SPINON)
                       CASE(0)
                          RHO = BO(L, K)
                       CASE(1)
                          RHO = RHOUP(L, K) + RHODOWN(L, K)
                          RHODIFF = (RHOUP(L, K) - RHODOWN(L, K))*(H2VECT(L) + H2VECT(K))
                       END SELECT

!                       RHO = X2HRHO(L, K)

                       IF (.NOT. PATH) THEN

                          ! Unroll loops and pre-compute

                          CALL DFDX_HS(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ALPHA, COSBETA, &
                               MYDFDA, MYDFDB, MYDFDR, &
                               SMYDFDA, SMYDFDB, SMYDFDR)

                          !
                          ! d/d_alpha
                          !
                          
                          FTMP_BOND(1) = FTMP_BOND(1) + RHO * &
                               (-RIJ(2) / MAGRP2 * MYDFDA)

                          FTMP_BOND(2) = FTMP_BOND(2) + RHO * &
                               (RIJ(1)/ MAGRP2 * MYDFDA)


                          FTMP_PULAY(1) = FTMP_PULAY(1) + X2HRHO(L, K) * &
                               (-RIJ(2) / MAGRP2 * SMYDFDA)
                          
                          FTMP_PULAY(2) = FTMP_PULAY(2) + X2HRHO(L, K) * &
                               (RIJ(1)/ MAGRP2 * SMYDFDA)


                          FTMP_COUL(1) = FTMP_COUL(1) + RHO * &
                               (-RIJ(2) / MAGRP2 * SMYDFDA)

                          FTMP_COUL(2) = FTMP_COUL(2) + RHO * &
                               (RIJ(1)/ MAGRP2 * SMYDFDA)

                          IF (SPINON .EQ. 1) THEN

                             FTMP_SPIN(1) = FTMP_SPIN(1) + RHODIFF * &
                                  (-RIJ(2) / MAGRP2 * SMYDFDA)
                             
                             FTMP_SPIN(2) = FTMP_SPIN(2) + RHODIFF * &
                                  (RIJ(1)/ MAGRP2 * SMYDFDA)

                          ENDIF


                          !
                          ! d/d_beta
                          !
                          
                          FTMP_BOND(1) = FTMP_BOND(1) + RHO * &
                               (((((RIJ(3) * RIJ(1)) / &
                               MAGR2)) / MAGRP) * MYDFDB)

                          FTMP_BOND(2) = FTMP_BOND(2) + RHO * &
                               (((((RIJ(3) * RIJ(2)) / &
                               MAGR2)) / MAGRP) * MYDFDB)

                          FTMP_BOND(3) = FTMP_BOND(3) - RHO * &
                               (((ONE - ((RIJ(3) * RIJ(3)) / &
                               MAGR2)) / MAGRP) * MYDFDB)


                          FTMP_PULAY(1) = FTMP_PULAY(1) + X2HRHO(L, K) * &
                               (((((RIJ(3) * RIJ(1)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)


                          FTMP_PULAY(2) = FTMP_PULAY(2) + X2HRHO(L, K) * &
                               (((((RIJ(3) * RIJ(2)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)
                               
                          FTMP_PULAY(3) = FTMP_PULAY(3) - X2HRHO(L, K) * &
                               (((ONE - ((RIJ(3) * RIJ(3)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)

                          
                          FTMP_COUL(1) = FTMP_COUL(1) + RHO * &
                               (((((RIJ(3) * RIJ(1)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)

                          FTMP_COUL(2) = FTMP_COUL(2) + RHO * &
                               (((((RIJ(3) * RIJ(2)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)

                          FTMP_COUL(3) = FTMP_COUL(3) - RHO * &
                               (((ONE - ((RIJ(3) * RIJ(3)) / &
                               MAGR2)) / MAGRP) * SMYDFDB)

                          IF (SPINON .EQ. 1) THEN
                             
                             FTMP_SPIN(1) = FTMP_SPIN(1) + RHODIFF * &
                                  (((((RIJ(3) * RIJ(1)) / &
                                  MAGR2)) / MAGRP) * SMYDFDB)
                             
                             FTMP_SPIN(2) = FTMP_SPIN(2) + RHODIFF * &
                                  (((((RIJ(3) * RIJ(2)) / &
                                  MAGR2)) / MAGRP) * SMYDFDB)
                             
                             FTMP_SPIN(3) = FTMP_SPIN(3) - RHODIFF * &
                                  (((ONE - ((RIJ(3) * RIJ(3)) / &
                                  MAGR2)) / MAGRP) * SMYDFDB)
                             
                          ENDIF


                          !
                          ! d/dR
                          !

                          FTMP_BOND(1) = FTMP_BOND(1) - RHO * DC(1) * &
                               MYDFDR

                          FTMP_BOND(2) = FTMP_BOND(2) - RHO * DC(2) * &
                               MYDFDR

                          FTMP_BOND(3) = FTMP_BOND(3) - RHO * DC(3) * &
                               MYDFDR


                          FTMP_PULAY(1) = FTMP_PULAY(1) - X2HRHO(L, K) * DC(1) * &
                               SMYDFDR

                          FTMP_PULAY(2) = FTMP_PULAY(2) - X2HRHO(L, K) * DC(2) * &
                               SMYDFDR

                          FTMP_PULAY(3) = FTMP_PULAY(3) - X2HRHO(L, K) * DC(3) * &
                               SMYDFDR
                          

                          FTMP_COUL(1) = FTMP_COUL(1) - RHO * DC(1) * &
                               SMYDFDR

                          FTMP_COUL(2) = FTMP_COUL(2) - RHO * DC(2) * &
                               SMYDFDR

                          FTMP_COUL(3) = FTMP_COUL(3) - RHO * DC(3) * &
                               SMYDFDR

                          IF (SPINON .EQ. 1) THEN

                             FTMP_SPIN(1) = FTMP_SPIN(1) - RHODIFF * DC(1) * &
                                  SMYDFDR

                             FTMP_SPIN(2) = FTMP_SPIN(2) - RHODIFF * DC(2) * &
                                  SMYDFDR
                             
                             FTMP_SPIN(3) = FTMP_SPIN(3) - RHODIFF * DC(3) * &
                                  SMYDFDR
                             
                          ENDIF
                             

                       ELSE

                          ! pathological configuration in which beta=0
                          ! or pi => alpha undefined

                          ! fixed: MJC 12/17/13

                          CALL DFDX_HS(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ZERO, COSBETA, &
                               MYDFDA, MYDFDB, MYDFDR, &
                               SMYDFDA, SMYDFDB, SMYDFDR)

                          MYDFDB = MYDFDB/MAGR
                          SMYDFDB = SMYDFDB/MAGR

                          
                          FTMP_BOND(1) = FTMP_BOND(1) - RHO * (COSBETA * MYDFDB)
                          
                          FTMP_BOND(3) = FTMP_BOND(3) - RHO * COSBETA * MYDFDR 

                          FTMP_PULAY(1) = FTMP_PULAY(1) - X2HRHO(L, K) * (COSBETA * SMYDFDB)

                          FTMP_COUL(1) = FTMP_COUL(1) - RHO * (COSBETA * SMYDFDB)

                          FTMP_PULAY(3) = FTMP_PULAY(3) - X2HRHO(L, K) * COSBETA * SMYDFDR

                          FTMP_COUL(3) = FTMP_COUL(3) - RHO * COSBETA * SMYDFDR

                          IF (SPINON .EQ. 1) THEN 
                             FTMP_SPIN(1) = FTMP_SPIN(1) - RHODIFF * (COSBETA * SMYDFDB)
                             FTMP_SPIN(3) = FTMP_SPIN(3) - RHODIFF * COSBETA * SMYDFDR
                          ENDIF

                          CALL DFDX_HS(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, PI/TWO, COSBETA, &
                               MYDFDA, MYDFDB, MYDFDR, &
                               SMYDFDA, SMYDFDB, SMYDFDR)
                          
                          MYDFDB = MYDFDB/MAGR
                          SMYDFDB = SMYDFDB/MAGR
                          
                          FTMP_BOND(2) = FTMP_BOND(2) - RHO * (COSBETA * MYDFDB)

                          FTMP_PULAY(2) = FTMP_PULAY(2) - X2HRHO(L, K) * (COSBETA * SMYDFDB)

                          FTMP_COUL(2) = FTMP_COUL(2) - RHO * (COSBETA * SMYDFDB)

                          IF (SPINON .EQ. 1) THEN
                             FTMP_SPIN(2) = FTMP_SPIN(2) - RHODIFF * (COSBETA * SMYDFDB)
                          ENDIF
                             
                       ENDIF

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           F(1,I) = F(1,I) + TWO*FTMP_BOND(1)
           F(2,I) = F(2,I) + TWO*FTMP_BOND(2)
           F(3,I) = F(3,I) + TWO*FTMP_BOND(3)

           VIRBOND(1) = VIRBOND(1) + RIJ(1) * FTMP_BOND(1)
           VIRBOND(2) = VIRBOND(2) + RIJ(2) * FTMP_BOND(2)
           VIRBOND(3) = VIRBOND(3) + RIJ(3) * FTMP_BOND(3)
           VIRBOND(4) = VIRBOND(4) + RIJ(1) * FTMP_BOND(2)
           VIRBOND(5) = VIRBOND(5) + RIJ(2) * FTMP_BOND(3)
           VIRBOND(6) = VIRBOND(6) + RIJ(3) * FTMP_BOND(1)

           FPUL(1,I) = FPUL(1,I) - TWO * FTMP_PULAY(1)
           FPUL(2,I) = FPUL(2,I) - TWO * FTMP_PULAY(2)
           FPUL(3,I) = FPUL(3,I) - TWO * FTMP_PULAY(3)

           VIRPUL(1) = VIRPUL(1) - RIJ(1) * FTMP_PULAY(1)
           VIRPUL(2) = VIRPUL(2) - RIJ(2) * FTMP_PULAY(2)
           VIRPUL(3) = VIRPUL(3) - RIJ(3) * FTMP_PULAY(3)
           VIRPUL(4) = VIRPUL(4) - RIJ(1) * FTMP_PULAY(2)
           VIRPUL(5) = VIRPUL(5) - RIJ(2) * FTMP_PULAY(3)
           VIRPUL(6) = VIRPUL(6) - RIJ(3) * FTMP_PULAY(1)

           IF (ELECTRO .EQ. 1) THEN

              FTMP_COUL = FTMP_COUL * ( HUBBARDU(ELEMPOINTER(J))*DELTAQ(J) + COULOMBV(J) &
                    +HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) + COULOMBV(I))

              FPUL(1,I) = FPUL(1,I) + FTMP_COUL(1)
              FPUL(2,I) = FPUL(2,I) + FTMP_COUL(2)
              FPUL(3,I) = FPUL(3,I) + FTMP_COUL(3)
              
              ! with the factor of 2...                  
              
              VIRPUL(1) = VIRPUL(1) + RIJ(1)*FTMP_COUL(1)/TWO
              VIRPUL(2) = VIRPUL(2) + RIJ(2)*FTMP_COUL(2)/TWO
              VIRPUL(3) = VIRPUL(3) + RIJ(3)*FTMP_COUL(3)/TWO
              VIRPUL(4) = VIRPUL(4) + RIJ(1)*FTMP_COUL(2)/TWO
              VIRPUL(5) = VIRPUL(5) + RIJ(2)*FTMP_COUL(3)/TWO
              VIRPUL(6) = VIRPUL(6) + RIJ(3)*FTMP_COUL(1)/TWO
              
           ELSE 

              FTMP_COUL = FTMP_COUL * (LCNSHIFT(I) + LCNSHIFT(J))

              FPUL(1,I) = FPUL(1,I) + FTMP_COUL(1)
              FPUL(2,I) = FPUL(2,I) + FTMP_COUL(2)
              FPUL(3,I) = FPUL(3,I) + FTMP_COUL(3)
              
              ! with the factor of 2...                  
              
              VIRPUL(1) = VIRPUL(1) + RIJ(1)*FTMP_COUL(1)/TWO
              VIRPUL(2) = VIRPUL(2) + RIJ(2)*FTMP_COUL(2)/TWO
              VIRPUL(3) = VIRPUL(3) + RIJ(3)*FTMP_COUL(3)/TWO
              VIRPUL(4) = VIRPUL(4) + RIJ(1)*FTMP_COUL(2)/TWO
              VIRPUL(5) = VIRPUL(5) + RIJ(2)*FTMP_COUL(3)/TWO
              VIRPUL(6) = VIRPUL(6) + RIJ(3)*FTMP_COUL(1)/TWO
              
           ENDIF

           IF (SPINON .EQ. 1) THEN
              
              FPUL(1,I) = FPUL(1,I) + FTMP_SPIN(1)
              FPUL(2,I) = FPUL(2,I) + FTMP_SPIN(2)
              FPUL(3,I) = FPUL(3,I) + FTMP_SPIN(3)
              
              ! with the factor of 2...                                           \
              
              VIRPUL(1) = VIRPUL(1) + RIJ(1)*FTMP_SPIN(1)/TWO
              VIRPUL(2) = VIRPUL(2) + RIJ(2)*FTMP_SPIN(2)/TWO
              VIRPUL(3) = VIRPUL(3) + RIJ(3)*FTMP_SPIN(3)/TWO
              VIRPUL(4) = VIRPUL(4) + RIJ(1)*FTMP_SPIN(2)/TWO
              VIRPUL(5) = VIRPUL(5) + RIJ(2)*FTMP_SPIN(3)/TWO
              VIRPUL(6) = VIRPUL(6) + RIJ(3)*FTMP_SPIN(1)/TWO
              

           ENDIF
           
        ENDIF

     ENDDO

  ENDDO

!$OMP END PARALLEL DO

  !  DO I= 1, NATS
  !     WRITE(6,10) I, FPUL(1,I), FPUL(2,I), FPUL(3,I)
  !  ENDDO

  !10 FORMAT(I4, 3F12.6)

  RETURN

END SUBROUTINE TBFORCES
