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

SUBROUTINE FSPINNONO_SP

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE NONOARRAY
  USE COULOMBARRAY
  USE NEBLISTARRAY
  USE SPINARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, KK, INDI, INDJ
  INTEGER :: NEWJ
  INTEGER :: PBCI, PBCJ, PBCK
  INTEGER :: SPININDI, SPININDJ 
  REAL(LATTEPREC) :: HSSS, HSPS, HPSS, HPPS, HPPP
  REAL(LATTEPREC) :: RIJ(3), DC(3)
  REAL(LATTEPREC) :: L, M, N, L2, M2, N2, LM, LN, MN, LMN
  REAL(LATTEPREC) :: DSSSDR(3), DSPSDR(3), DPSSDR(3), DPPSDR(3), DPPPDR(3)
  REAL(LATTEPREC) :: MAGR, INVR, FTMPS(3), FTMPP(3), FTMP(3)
  REAL(LATTEPREC) :: PPSMPPP, PPSUBINVR
  REAL(LATTEPREC) :: VIRTMP(6), VIRTMPS(6), VIRTMPP(6)
  CHARACTER(LEN=2) :: BASISI, BASISJ
  IF (EXISTERROR) RETURN

  !  INDI = 0

  FSSPIN = ZERO

  VIRSSPIN = ZERO

  !
  ! We are computing the contribution to the forces from the 
  ! S dependence of the Mulliken spin densities when we use 
  ! the non-orthogonal basis

  DO I = 1, NATS

     BASISI = BASIS(ELEMPOINTER(I))
     INDI = MATINDLIST(I)
     SPININDI = SPININDLIST(I)
     ! Loop over all neighbors of I

     DO NEWJ = 1, TOTNEBTB(I)

        J = NEBTB(1, NEWJ, I)
        PBCI = NEBTB(2, NEWJ, I)
        PBCJ = NEBTB(3, NEWJ, I)
        PBCK = NEBTB(4, NEWJ, I)        

        INDJ = MATINDLIST(J)
        SPININDJ = SPININDLIST(J)

        BASISJ = BASIS(ELEMPOINTER(J))

        RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1) - CR(1,I)

        RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2) - CR(2,I)

        RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3) - CR(3,I)

        MAGR = SQRT(RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3))

        INVR = ONE/MAGR

        FTMPS = ZERO
        FTMPP = ZERO

        ! 
        ! Direction cosines (DC)
        !

        DC = RIJ/MAGR

        L = DC(1)
        M = DC(2)
        N = DC(3)

        ! Let's compute rho * dS/dR

        IF (BASISI .EQ. "s") THEN

           IF (BASISJ .EQ. "s") THEN

              DO K = 1, NOINT
                 IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                      ATELE(J) .EQ. ELE2(K)) .OR. &
                      (ATELE(I) .EQ. ELE2(K) .AND. &
                      ATELE(J) .EQ. ELE1(K))) THEN

                    IF (BTYPE(K) .EQ. "sss") THEN

                       CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                    ENDIF

                 ENDIF
              ENDDO

              FTMPS = FTMPS -  DSSSDR*(RHOUP(INDI+1, INDJ+1) - &
                   RHODOWN(INDI+1, INDJ+1))

           ELSEIF (BASISJ .EQ. "sp") THEN

              DO K = 1, NOINT

                 IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                      ATELE(J) .EQ. ELE2(K)) .OR. &
                      (ATELE(I) .EQ. ELE2(K) .AND. &
                      ATELE(J) .EQ. ELE1(K))) THEN

                    IF (BTYPE(K) .EQ. "sss") THEN

                       CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                    ELSEIF (BTYPE(K) .EQ. "sps") THEN

                       CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSPS, DSPSDR)

                    ENDIF
                 ENDIF
              ENDDO

              L2 = L*L
              M2 = M*M
              N2 = N*N
              LM = L*M
              LN = L*N
              MN = M*N

              ! E_s1,s2

              FTMPS = FTMPS - DSSSDR*(RHOUP(INDI+1, INDJ+1) - &
                   RHODOWN(INDI+1, INDJ+1))

              ! E_s1,x2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(1) + (L2 - ONE)*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(2) + LM*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(3) + LN*INVR*HSPS)

              ! E_s1,y2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(1) + LM*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(2) + (M2 - ONE)*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(3) + MN*INVR*HSPS)

              ! E_s1,z2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(1) + LN*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(2) + MN*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(3) + (N2 - ONE)*INVR*HSPS)

           ENDIF

        ELSEIF (BASISI .EQ. "sp") THEN

           IF (BASISJ .EQ. "s") THEN

              DO K = 1, NOINT

                 IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                      ATELE(J) .EQ. ELE2(K)) .OR. &
                      (ATELE(I) .EQ. ELE2(K) .AND. &
                      ATELE(J) .EQ. ELE1(K))) THEN

                    IF (BTYPE(K) .EQ. "sss") THEN

                       CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                    ELSEIF (BTYPE(K) .EQ. "sps") THEN

                       CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPSS, DPSSDR)

                       HPSS = -HPSS
                       DPSSDR = -DPSSDR

                    ENDIF
                 ENDIF
              ENDDO

              L2 = L*L
              M2 = M*M
              N2 = N*N
              LM = L*M
              LN = L*N
              MN = M*N

              ! E_s1,s2

              FTMPS = FTMPS - DSSSDR*(RHOUP(INDI+1, INDJ+1) - &
                   RHODOWN(INDI+1, INDJ+1))

              ! E_x1,s2

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1)) * &
                   (L*DPSSDR(1) + (L2 - ONE)*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1)) * &
                   (L*DPSSDR(2) + LM*INVR*HPSS)

              FTMPS(3) = FTMPS(3) - (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1)) * &
                   (L*DPSSDR(3) + LN*INVR*HPSS)

              ! E_y1,s2

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(1) + LM*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(2) + (M2 - ONE)*INVR*HPSS)

              FTMPS(3) = FTMPS(3) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(3) + MN*INVR*HPSS)

              ! E_z1,s2

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(1) + LN*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(2) + MN*INVR*HPSS)

              FTMPS(3) = FTMPS(3) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(3) + (N2 - ONE)*INVR*HPSS)

           ELSEIF (BASISJ .EQ. "sp") THEN

              IF (ATELE(I) .EQ. ATELE(J)) THEN

                 DO K = 1, NOINT

                    IF (ATELE(I) .EQ. ELE1(K) .AND. &
                         ATELE(J) .EQ. ELE2(K)) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN  

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                       ELSEIF (BTYPE(K) .EQ. "sps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSPS, DSPSDR)

                          DPSSDR = -DSPSDR

                          HPSS = -HSPS

                       ELSEIF (BTYPE(K) .EQ. "pps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPS, DPPSDR)

                       ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPP, DPPPDR)

                       ENDIF
                    ENDIF
                 ENDDO

              ELSEIF (ATELE(I) .NE. ATELE(J)) THEN

                 DO K = 1, NOINT

                    IF (ATELE(I) .EQ. ELE1(K) .AND. &
                         ATELE(J) .EQ. ELE2(K)) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN                 

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                       ELSEIF (BTYPE(K) .EQ. "sps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSPS, DSPSDR)

                       ELSEIF (BTYPE(K) .EQ. "pps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPS, DPPSDR)

                       ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPP, DPPPDR)

                       ENDIF

                    ELSEIF (ATELE(I) .EQ. ELE2(K) .AND. &
                         ATELE(J) .EQ. ELE1(K)) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HSSS, DSSSDR)

                       ELSEIF (BTYPE(K) .EQ. "sps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPSS, DPSSDR)

                          DPSSDR = -DPSSDR
                          HPSS = -HPSS

                       ELSEIF (BTYPE(K) .EQ. "pps") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPS, DPPSDR)

                       ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                          CALL DUNIVSCALE_SUB(MAGR, OVERL(:,K), DC, HPPP, DPPPDR)

                       ENDIF

                    ENDIF
                 ENDDO

              ENDIF

              PPSMPPP = HPPS - HPPP
              PPSUBINVR = PPSMPPP * INVR

              L2 = L*L
              M2 = M*M
              N2 = N*N
              LM = L*M
              LN = L*N
              MN = M*N
              LMN = LM*N

              ! E_s1,s2

              FTMPS = FTMPS - DSSSDR*(RHOUP(INDI+1, INDJ+1) - &
                   RHODOWN(INDI+1, INDJ+1))

              ! E_s1,x2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(1) + (L2 - ONE)*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(2) + LM*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+2) - &
                   RHODOWN(INDI+1, INDJ+2)) * &
                   (L*DSPSDR(3) + LN*INVR*HSPS)

              ! E_s1,y2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(1) + LM*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(2) + (M2 - ONE)*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+3) - &
                   RHODOWN(INDI+1, INDJ+3)) * &
                   (M*DSPSDR(3) + MN*INVR*HSPS)

              ! E_s1,z2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(1) + LN*INVR*HSPS)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(2) + MN*INVR*HSPS)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+1, INDJ+4) - &
                   RHODOWN(INDI+1, INDJ+4)) * &
                   (N*DSPSDR(3) + (N2 - ONE)*INVR*HSPS)

              ! E_x1,s2  

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1)) * &
                   (L*DPSSDR(1) + (L2 - ONE)*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1)) * &
                   (L*DPSSDR(2) + LM*INVR*HPSS)

              FTMPS(3) = FTMPS(3) -  (RHOUP(INDI+2, INDJ+1) - &
                   RHODOWN(INDI+2, INDJ+1))* &
                   (L*DPSSDR(3) + LN*INVR*HPSS)

              ! E_x1,x2

              FTMPP(1) = FTMPP(1) -  (RHOUP(INDI+2, INDJ+2) - &
                   RHODOWN(INDI+2, INDJ+2)) * &
                   (L2*DPPSDR(1) + (ONE - L2)*DPPPDR(1) + &
                   TWO*L*(L2 - ONE)*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+2, INDJ+2) - &
                   RHODOWN(INDI+2, INDJ+2)) * &
                   (L2*DPPSDR(2) + (ONE - L2)*DPPPDR(2) + &
                   TWO*L2*M*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+2, INDJ+2) - &
                   RHODOWN(INDI+2, INDJ+2)) * &
                   (L2*DPPSDR(3) + (ONE - L2)*DPPPDR(3) + &
                   TWO*L2*N*PPSUBINVR)

              ! E_x1,y2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+2, INDJ+3) - &
                   RHODOWN(INDI+2, INDJ+3)) * &
                   (LM*(DPPSDR(1) - DPPPDR(1)) + &
                   M*(TWO*L2 - ONE)*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+2, INDJ+3) - &
                   RHODOWN(INDI+2, INDJ+3)) * &
                   (LM*(DPPSDR(2) - DPPPDR(2)) + &
                   L*(TWO*M2 - ONE)*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+2, INDJ+3) - &
                   RHODOWN(INDI+2, INDJ+3)) * &
                   (LM*(DPPSDR(3) - DPPPDR(3)) + &
                   TWO*LMN*PPSUBINVR)

              ! E_x1,z2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+2, INDJ+4) - &
                   RHODOWN(INDI+2, INDJ+4)) * &
                   (LN*(DPPSDR(1) - DPPPDR(1)) + &
                   N*(TWO*L2 - ONE)*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+2, INDJ+4) - &
                   RHODOWN(INDI+2, INDJ+4))* &
                   (LN*(DPPSDR(2) - DPPPDR(2)) + &
                   TWO*LMN*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+2, INDJ+4) - &
                   RHODOWN(INDI+2, INDJ+4)) * &
                   (LN*(DPPSDR(3) - DPPPDR(3)) + &
                   L*(TWO*N2 - ONE)*PPSUBINVR)

              ! E_y1,s2

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(1) + LM*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(2) + (M2 - ONE)*INVR*HPSS)

              FTMPS(3) = FTMPS(3) - (RHOUP(INDI+3, INDJ+1) - &
                   RHODOWN(INDI+3, INDJ+1)) * &
                   (M*DPSSDR(3) + MN*INVR*HPSS)

              ! E_y1,x2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+3, INDJ+2) - &
                   RHODOWN(INDI+3, INDJ+2)) * &
                   (LM*(DPPSDR(1) - DPPPDR(1)) + &
                   M*(TWO*L2 - ONE)*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+3, INDJ+2) - &
                   RHODOWN(INDI+3, INDJ+2)) * &
                   (LM*(DPPSDR(2) - DPPPDR(2)) + &
                   L*(TWO*M2 - ONE)*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+3, INDJ+2) - &
                   RHODOWN(INDI+3, INDJ+2)) * &
                   (LM*(DPPSDR(3) - DPPPDR(3)) + &
                   TWO*LMN*PPSUBINVR)

              ! E_y1,y2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+3, INDJ+3) - &
                   RHODOWN(INDI+3, INDJ+3)) * &
                   (M2*DPPSDR(1) + (ONE - M2)*DPPPDR(1) + &
                   TWO*L*M2*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+3, INDJ+3) - &
                   RHODOWN(INDI+3, INDJ+3)) * &
                   (M2*DPPSDR(2) + (ONE - M2)*DPPPDR(2) + &
                   TWO*M*(M2 - ONE)*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+3, INDJ+3) - &
                   RHODOWN(INDI+3, INDJ+3)) * &
                   (M2*DPPSDR(3) + (ONE - M2)*DPPPDR(3) + &
                   TWO*N*M2*PPSUBINVR)

              ! E_y1,z2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+3, INDJ+4) - &
                   RHODOWN(INDI+3, INDJ+4)) * &
                   (MN*(DPPSDR(1) - DPPPDR(1)) + &
                   TWO*LMN*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+3, INDJ+4) - &
                   RHODOWN(INDI+3, INDJ+4)) * &
                   (MN*(DPPSDR(2) - DPPPDR(2)) + &
                   N*(TWO*M2 - ONE)*PPSUBINVR)

              FTMPP(3) = FTMPP(3) -  (RHOUP(INDI+3, INDJ+4) - &
                   RHODOWN(INDI+3, INDJ+4))* &
                   (MN*(DPPSDR(3) - DPPPDR(3)) + &
                   M*(TWO*N2 - ONE)*PPSUBINVR)

              ! E_z1,s2

              FTMPS(1) = FTMPS(1) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(1) + LN*INVR*HPSS)

              FTMPS(2) = FTMPS(2) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(2) + MN*INVR*HPSS)

              FTMPS(3) = FTMPS(3) - (RHOUP(INDI+4, INDJ+1) - &
                   RHODOWN(INDI+4, INDJ+1)) * &
                   (N*DPSSDR(3) + (N2 - ONE)*INVR*HPSS)

              ! E_z1,x2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+4, INDJ+2) - &
                   RHODOWN(INDI+4, INDJ+2)) * &
                   (LN*(DPPSDR(1) - DPPPDR(1)) + &
                   N*(TWO*L2 - ONE)*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+4, INDJ+2) - &
                   RHODOWN(INDI+4, INDJ+2)) * &
                   (LN*(DPPSDR(2) - DPPPDR(2)) + &
                   TWO*LMN*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+4, INDJ+2) - &
                   RHODOWN(INDI+4, INDJ+2)) * &
                   (LN*(DPPSDR(3) - DPPPDR(3)) + &
                   L*(TWO*N2 - ONE)*PPSUBINVR)

              ! E_z1,y2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+4, INDJ+3) - &
                   RHODOWN(INDI+4, INDJ+3)) * &
                   (MN*(DPPSDR(1) - DPPPDR(1)) + &
                   TWO*LMN*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+4, INDJ+3) - &
                   RHODOWN(INDI+4, INDJ+3)) * &
                   (MN*(DPPSDR(2) - DPPPDR(2)) + &
                   N*(TWO*M2 - ONE)*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+4, INDJ+3) - &
                   RHODOWN(INDI+4, INDJ+3)) * &
                   (MN*(DPPSDR(3) - DPPPDR(3)) + &
                   M*(TWO*N2 - ONE)*PPSUBINVR)

              ! E_z1,z2

              FTMPP(1) = FTMPP(1) - (RHOUP(INDI+4, INDJ+4) - &
                   RHODOWN(INDI+4, INDJ+4)) * &
                   (N2*DPPSDR(1) + (ONE - N2)*DPPPDR(1) + &
                   TWO*L*N2*PPSUBINVR)

              FTMPP(2) = FTMPP(2) - (RHOUP(INDI+4, INDJ+4) - &
                   RHODOWN(INDI+4, INDJ+4)) * &
                   (N2*DPPSDR(2) + (ONE - N2)*DPPPDR(2) + &
                   TWO*M*N2*PPSUBINVR)

              FTMPP(3) = FTMPP(3) - (RHOUP(INDI+4, INDJ+4) - &
                   RHODOWN(INDI+4, INDJ+4)) * &
                   (N2*DPPSDR(3) + (ONE - N2)*DPPPDR(3) + &
                   TWO*N*(N2 - ONE)*PPSUBINVR)

           ENDIF

        ENDIF

        IF (BASISJ .EQ. "s") THEN

           FTMP = FTMPS*DELTASPIN(SPININDJ+1)*WSS(ELEMPOINTER(J))

        ELSEIF (BASISJ .EQ. "sp") THEN

           FTMP = FTMPS*DELTASPIN(SPININDJ+1)*WSS(ELEMPOINTER(J)) + &
                FTMPP*DELTASPIN(SPININDJ+2)*WPP(ELEMPOINTER(J))

        ENDIF

        IF (BASISI .EQ. "s") THEN

           FTMP = FTMP + FTMPS*DELTASPIN(SPININDI+1)*WSS(ELEMPOINTER(I))

        ELSEIF (BASISI .EQ. "sp") THEN

           FTMP = FTMP + FTMPS*DELTASPIN(SPININDI+1)*WSS(ELEMPOINTER(I)) + &
                FTMPP*DELTASPIN(SPININDI+2)*WPP(ELEMPOINTER(I))

        ENDIF

        FSSPIN(1,I) = FSSPIN(1,I) + FTMP(1)
        FSSPIN(2,I) = FSSPIN(2,I) + FTMP(2)
        FSSPIN(3,I) = FSSPIN(3,I) + FTMP(3)

        ! with the factor of 2...

        VIRSSPIN(1) = VIRSSPIN(1) + RIJ(1)*FTMP(1)
        VIRSSPIN(2) = VIRSSPIN(2) + RIJ(2)*FTMP(2)
        VIRSSPIN(3) = VIRSSPIN(3) + RIJ(3)*FTMP(3)
        VIRSSPIN(4) = VIRSSPIN(4) + RIJ(1)*FTMP(2)
        VIRSSPIN(5) = VIRSSPIN(5) + RIJ(2)*FTMP(3)
        VIRSSPIN(6) = VIRSSPIN(6) + RIJ(3)*FTMP(1)

     ENDDO

  ENDDO

  VIRSSPIN = HALF*VIRSSPIN

  RETURN

END SUBROUTINE FSPINNONO_SP
