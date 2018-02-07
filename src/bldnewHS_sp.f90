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

SUBROUTINE BLDNEWHS_SP

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NEBLISTARRAY
  USE XBOARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE UNIVARRAY
  USE MYPRECISION

#ifdef PROGRESSON
  USE GENXPROGRESS
#endif

  IMPLICIT NONE

  INTEGER :: I, J, NEWJ, K, KK, SUBI, NNZ
  INTEGER :: MYINDEX, INDI, INDJ
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: HPPS, HPPP, PPSMPPP
  REAL(LATTEPREC) :: L, M, N, RIJ(3), MAGR, MAGR2, RCUTTB

  !
  ! The sp H-matrix elements so we don't confuse ourselves...
  !

  REAL(LATTEPREC) :: SISJ = ZERO 
  REAL(LATTEPREC) :: SIPXJ = ZERO, SIPYJ = ZERO, SIPZJ = ZERO
  REAL(LATTEPREC) :: PXISJ = ZERO, PYISJ = ZERO, PZISJ = ZERO
  REAL(LATTEPREC) :: PXPX = ZERO, PXPY = ZERO, PXPZ = ZERO
  REAL(LATTEPREC) :: PYPX = ZERO, PYPY = ZERO, PYPZ = ZERO
  REAL(LATTEPREC) :: PZPX = ZERO, PZPY = ZERO, PZPZ = ZERO

  REAL(LATTEPREC) :: PISJ, OLPXISJ, OLPYISJ, OLPZISJ
  REAL(LATTEPREC) :: OLSIPXJ, OLSIPYJ, OLSIPZJ, OLSISJ, SIPJ
  REAL(LATTEPREC) :: SPPP, SPPS 

  !

  CHARACTER(LEN=2) :: BASISI, BASISJ
  IF (EXISTERROR) RETURN

  H = ZERO

  MYINDEX = 0

  DO I = 1, NATS

     IF (BASIS(ELEMPOINTER(I)) .EQ. "sp") THEN

        MYINDEX = MYINDEX + 1
        H(MYINDEX, MYINDEX) = HES(ELEMPOINTER(I))
        MYINDEX = MYINDEX + 1
        H(MYINDEX, MYINDEX) = HEP(ELEMPOINTER(I))
        MYINDEX = MYINDEX + 1
        H(MYINDEX, MYINDEX) = HEP(ELEMPOINTER(I))
        MYINDEX = MYINDEX + 1
        H(MYINDEX, MYINDEX) = HEP(ELEMPOINTER(I))

     ELSEIF (BASIS(ELEMPOINTER(I)) .EQ. "s") THEN

        MYINDEX = MYINDEX + 1

        H(MYINDEX, MYINDEX) = HES(ELEMPOINTER(I))

     ENDIF

  ENDDO

  IF (BASISTYPE .EQ. "NONORTHO") THEN

     SMAT = ZERO

     DO I = 1, HDIM
        SMAT(I,I) = ONE
     ENDDO

  ENDIF

  !
  ! In the following, just to be clear, 
  !
  ! SSS = S S SIGMA
  ! SPS = S P SIGMA
  ! PSS = P S SIGMA
  ! PPS = P P SIGMA
  ! PPP = P P PI,
  ! 
  ! are the fundamental Slater-Koster bond integrals
  !

  !$OMP PARALLEL DO DEFAULT (NONE) &
  !$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &    
  !$OMP SHARED(CR, BOX, H, SMAT, NOINT, ATELE, ELE1, ELE2) &           
  !$OMP SHARED(BOND, OVERL, MATINDLIST, BTYPE, BASISTYPE) &
  !$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &
  !$OMP PRIVATE(RIJ, MAGR, MAGR2, SISJ, SIPXJ, SIPYJ, SIPZJ, PXISJ, PYISJ, PZISJ) &  
  !$OMP PRIVATE(PXPX, PXPY, PXPZ, PYPX, PYPY, PYPZ, PZPX, PZPY, PZPZ) &
  !$OMP PRIVATE(PISJ, OLPXISJ, OLPYISJ, OLPZISJ, OLSIPXJ, OLSIPYJ, OLSIPZJ, OLSISJ)&
  !$OMP PRIVATE(SIPJ, SPPP, SPPS, L, M, N, HPPS, HPPP, PPSMPPP)


  DO I = 1, NATS

     BASISI = BASIS(ELEMPOINTER(I))
     INDI = MATINDLIST(I)

     DO NEWJ = 1, TOTNEBTB(I)

        !
        ! Getting neighbors from the TB neighborlist
        !

        J = NEBTB(1, NEWJ, I)

        IF (J .GE. I) THEN !bug fix by MJC

           PBCI = NEBTB(2, NEWJ, I)
           PBCJ = NEBTB(3, NEWJ, I)
           PBCK = NEBTB(4, NEWJ, I)

           INDJ = MATINDLIST(J)

           !
           !     Which orbitals does J have? (s, sp, etc.)?
           !

           BASISJ = BASIS(ELEMPOINTER(J))

           RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
                REAL(PBCK)*BOX(3,1) - CR(1,I)

           RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
                REAL(PBCK)*BOX(3,2) - CR(2,I)

           RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
                REAL(PBCK)*BOX(3,3) - CR(3,I)

           MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

           MAGR = SQRT(MAGR2)

           !
           ! Direction cosines
           !

           L = RIJ(1)/MAGR
           M = RIJ(2)/MAGR
           N = RIJ(3)/MAGR

           !
           ! 5/3/10: Major bug fix by Ed Sanville for when 
           ! the periodic cell measures less than the cut-off
           ! for the bond integrals.
           !

           IF (BASISI .EQ. "s") THEN

              IF (BASISJ .EQ. "s") THEN

                 ! SISJ 

                 DO K = 1, NOINT
                    IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                         ATELE(J) .EQ. ELE2(K)) .OR. &
                         (ATELE(I) .EQ. ELE2(K) .AND. &
                         ATELE(J) .EQ. ELE1(K))) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN

                          CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                          IF (BASISTYPE .EQ. "NONORTHO") &
                               CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                       ENDIF

                    ENDIF
                 ENDDO

                 H(INDI+1, INDJ+1) = H(INDI+1, INDJ+1) + SISJ

                 IF (BASISTYPE .EQ. "NONORTHO") &
                      SMAT(INDI+1, INDJ+1) = SMAT(INDI+1, INDJ+1) + OLSISJ

              ELSEIF (BASISJ .EQ. "sp") THEN

                 ! SISJ, SIPXJ, SIPYJ, SIPZJ

                 DO K = 1, NOINT

                    IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                         ATELE(J) .EQ. ELE2(K)) .OR. &
                         (ATELE(I) .EQ. ELE2(K) .AND. &
                         ATELE(J) .EQ. ELE1(K))) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN

                          CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                          IF (BASISTYPE .EQ. "NONORTHO") &
                               CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                       ELSEIF (BTYPE(K) .EQ. "sps") THEN

                          CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SIPJ)

                          SIPXJ = L * SIPJ
                          SIPYJ = M * SIPJ
                          SIPZJ = N * SIPJ

                          IF (BASISTYPE .EQ. "NONORTHO") THEN

                             CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SIPJ)

                             OLSIPXJ = L * SIPJ
                             OLSIPYJ = M * SIPJ
                             OLSIPZJ = N * SIPJ

                          ENDIF

                       ENDIF
                    ENDIF
                 ENDDO

                 H(INDI+1, INDJ+1) = H(INDI+1, INDJ+1) + SISJ
                 H(INDI+1, INDJ+2) = H(INDI+1, INDJ+2) + SIPXJ
                 H(INDI+1, INDJ+3) = H(INDI+1, INDJ+3) + SIPYJ
                 H(INDI+1, INDJ+4) = H(INDI+1, INDJ+4) + SIPZJ

                 IF (BASISTYPE .EQ. "NONORTHO") THEN

                    SMAT(INDI+1, INDJ+1) = SMAT(INDI+1, INDJ+1) + OLSISJ
                    SMAT(INDI+1, INDJ+2) = SMAT(INDI+1, INDJ+2) + OLSIPXJ
                    SMAT(INDI+1, INDJ+3) = SMAT(INDI+1, INDJ+3) + OLSIPYJ
                    SMAT(INDI+1, INDJ+4) = SMAT(INDI+1, INDJ+4) + OLSIPZJ

                 ENDIF

              ENDIF

           ELSEIF (BASISI .EQ. "sp") THEN

              IF (BASISJ .EQ. "s") THEN

                 ! SISJ, PXISJ, PYISJ, PZISJ

                 DO K = 1, NOINT

                    IF ((ATELE(I) .EQ. ELE1(K) .AND. &
                         ATELE(J) .EQ. ELE2(K)) .OR. &
                         (ATELE(I) .EQ. ELE2(K) .AND. &
                         ATELE(J) .EQ. ELE1(K))) THEN

                       IF (BTYPE(K) .EQ. "sss") THEN

                          CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                          IF (BASISTYPE .EQ. "NONORTHO") &
                               CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                       ELSEIF (BTYPE(K) .EQ. "sps") THEN

                          CALL UNIVSCALE_SUB(MAGR, BOND(:,K), PISJ)

                          PXISJ = -L * PISJ
                          PYISJ = -M * PISJ
                          PZISJ = -N * PISJ

                          IF (BASISTYPE .EQ. "NONORTHO") THEN

                             CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), PISJ)

                             OLPXISJ = -L * PISJ
                             OLPYISJ = -M * PISJ
                             OLPZISJ = -N * PISJ

                          ENDIF

                       ENDIF
                    ENDIF
                 ENDDO

                 H(INDI+1, INDJ+1) = H(INDI+1, INDJ+1) + SISJ
                 H(INDI+2, INDJ+1) = H(INDI+2, INDJ+1) + PXISJ
                 H(INDI+3, INDJ+1) = H(INDI+3, INDJ+1) + PYISJ
                 H(INDI+4, INDJ+1) = H(INDI+4, INDJ+1) + PZISJ

                 IF (BASISTYPE .EQ. "NONORTHO") THEN

                    SMAT(INDI+1, INDJ+1) = SMAT(INDI+1, INDJ+1) + OLSISJ
                    SMAT(INDI+2, INDJ+1) = SMAT(INDI+2, INDJ+1) + OLPXISJ
                    SMAT(INDI+3, INDJ+1) = SMAT(INDI+3, INDJ+1) + OLPYISJ
                    SMAT(INDI+4, INDJ+1) = SMAT(INDI+4, INDJ+1) + OLPZISJ

                 ENDIF

              ELSEIF (BASISJ .EQ. "sp") THEN

                 ! element I = element J is a bit simpler 

                 IF (ATELE(I) .EQ. ATELE(J)) THEN

                    DO K = 1, NOINT

                       IF (ATELE(I) .EQ. ELE1(K) .AND. &
                            ATELE(J) .EQ. ELE2(K)) THEN

                          IF (BTYPE(K) .EQ. "sss") THEN                 

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                          ELSEIF (BTYPE(K) .EQ. "sps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SIPJ)

                             SIPXJ = L * SIPJ
                             SIPYJ = M * SIPJ
                             SIPZJ = N * SIPJ

                             PXISJ = -SIPXJ
                             PYISJ = -SIPYJ
                             PZISJ = -SIPZJ

                             IF (BASISTYPE .EQ. "NONORTHO") THEN

                                CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SIPJ)

                                OLSIPXJ = L * SIPJ
                                OLSIPYJ = M * SIPJ
                                OLSIPZJ = N * SIPJ

                                OLPXISJ = -OLSIPXJ
                                OLPYISJ = -OLSIPYJ
                                OLPZISJ = -OLSIPZJ

                             ENDIF

                          ELSEIF (BTYPE(K) .EQ. "pps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPS)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPS)

                          ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPP)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPP)

                          ENDIF
                       ENDIF
                    ENDDO

                 ELSEIF (ATELE(I) .NE. ATELE(J)) THEN

                    DO K = 1, NOINT

                       IF (ATELE(I) .EQ. ELE1(K) .AND. &
                            ATELE(J) .EQ. ELE2(K)) THEN

                          IF (BTYPE(K) .EQ. "sss") THEN                 

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                          ELSEIF (BTYPE(K) .EQ. "sps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SIPJ)

                             SIPXJ = L * SIPJ
                             SIPYJ = M * SIPJ
                             SIPZJ = N * SIPJ

                             IF (BASISTYPE .EQ. "NONORTHO") THEN

                                CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SIPJ)

                                OLSIPXJ = L * SIPJ
                                OLSIPYJ = M * SIPJ
                                OLSIPZJ = N * SIPJ

                             ENDIF

                          ELSEIF (BTYPE(K) .EQ. "pps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPS)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPS)

                          ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPP)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPP)

                          ENDIF

                       ELSEIF (ATELE(I) .EQ. ELE2(K) .AND. &
                            ATELE(J) .EQ. ELE1(K)) THEN

                          IF (BTYPE(K) .EQ. "sss") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SISJ)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), OLSISJ)

                          ELSEIF (BTYPE(K) .EQ. "sps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), SIPJ)

                             PXISJ = -L * SIPJ
                             PYISJ = -M * SIPJ
                             PZISJ = -N * SIPJ

                             IF (BASISTYPE .EQ. "NONORTHO") THEN

                                CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SIPJ)

                                OLPXISJ = -L * SIPJ
                                OLPYISJ = -M * SIPJ
                                OLPZISJ = -N * SIPJ

                             ENDIF

                          ELSEIF (BTYPE(K) .EQ. "pps") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPS)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPS)

                          ELSEIF (BTYPE(K) .EQ. "ppp") THEN

                             CALL UNIVSCALE_SUB(MAGR, BOND(:,K), HPPP)

                             IF (BASISTYPE .EQ. "NONORTHO") &
                                  CALL UNIVSCALE_SUB(MAGR, OVERL(:,K), SPPP)

                          ENDIF

                       ENDIF
                    ENDDO

                 ENDIF

                 PPSMPPP = HPPS - HPPP

                 PXPX = HPPP + L*L*PPSMPPP
                 PXPY = L*M*PPSMPPP
                 PXPZ = L*N*PPSMPPP
                 PYPX = M*L*PPSMPPP
                 PYPY = HPPP + M*M*PPSMPPP
                 PYPZ = M*N*PPSMPPP
                 PZPX = N*L*PPSMPPP
                 PZPY = N*M*PPSMPPP
                 PZPZ = HPPP + N*N*PPSMPPP

                 H(INDI+1,INDJ+1) = H(INDI+1,INDJ+1) + SISJ
                 H(INDI+1,INDJ+2) = H(INDI+1,INDJ+2) + SIPXJ
                 H(INDI+1,INDJ+3) = H(INDI+1,INDJ+3) + SIPYJ
                 H(INDI+1,INDJ+4) = H(INDI+1,INDJ+4) + SIPZJ

                 H(INDI+2,INDJ+1) = H(INDI+2,INDJ+1) + PXISJ
                 H(INDI+2,INDJ+2) = H(INDI+2,INDJ+2) + PXPX
                 H(INDI+2,INDJ+3) = H(INDI+2,INDJ+3) + PXPY
                 H(INDI+2,INDJ+4) = H(INDI+2,INDJ+4) + PXPZ
                 H(INDI+3,INDJ+1) = H(INDI+3,INDJ+1) + PYISJ
                 H(INDI+3,INDJ+2) = H(INDI+3,INDJ+2) + PYPX
                 H(INDI+3,INDJ+3) = H(INDI+3,INDJ+3) + PYPY
                 H(INDI+3,INDJ+4) = H(INDI+3,INDJ+4) + PYPZ
                 H(INDI+4,INDJ+1) = H(INDI+4,INDJ+1) + PZISJ
                 H(INDI+4,INDJ+2) = H(INDI+4,INDJ+2) + PZPX
                 H(INDI+4,INDJ+3) = H(INDI+4,INDJ+3) + PZPY
                 H(INDI+4,INDJ+4) = H(INDI+4,INDJ+4) + PZPZ

                 IF (BASISTYPE .EQ. "NONORTHO") THEN

                    PPSMPPP = SPPS - SPPP

                    PXPX = SPPP + L*L*PPSMPPP
                    PXPY = L*M*PPSMPPP
                    PXPZ = L*N*PPSMPPP
                    PYPX = M*L*PPSMPPP
                    PYPY = SPPP + M*M*PPSMPPP
                    PYPZ = M*N*PPSMPPP
                    PZPX = N*L*PPSMPPP
                    PZPY = N*M*PPSMPPP
                    PZPZ = SPPP + N*N*PPSMPPP

                    SMAT(INDI+1,INDJ+1) = SMAT(INDI+1,INDJ+1) + OLSISJ
                    SMAT(INDI+1,INDJ+2) = SMAT(INDI+1,INDJ+2) + OLSIPXJ
                    SMAT(INDI+1,INDJ+3) = SMAT(INDI+1,INDJ+3) + OLSIPYJ
                    SMAT(INDI+1,INDJ+4) = SMAT(INDI+1,INDJ+4) + OLSIPZJ

                    SMAT(INDI+2,INDJ+1) = SMAT(INDI+2,INDJ+1) + OLPXISJ
                    SMAT(INDI+2,INDJ+2) = SMAT(INDI+2,INDJ+2) + PXPX
                    SMAT(INDI+2,INDJ+3) = SMAT(INDI+2,INDJ+3) + PXPY
                    SMAT(INDI+2,INDJ+4) = SMAT(INDI+2,INDJ+4) + PXPZ
                    SMAT(INDI+3,INDJ+1) = SMAT(INDI+3,INDJ+1) + OLPYISJ
                    SMAT(INDI+3,INDJ+2) = SMAT(INDI+3,INDJ+2) + PYPX
                    SMAT(INDI+3,INDJ+3) = SMAT(INDI+3,INDJ+3) + PYPY
                    SMAT(INDI+3,INDJ+4) = SMAT(INDI+3,INDJ+4) + PYPZ
                    SMAT(INDI+4,INDJ+1) = SMAT(INDI+4,INDJ+1) + OLPZISJ
                    SMAT(INDI+4,INDJ+2) = SMAT(INDI+4,INDJ+2) + PZPX
                    SMAT(INDI+4,INDJ+3) = SMAT(INDI+4,INDJ+3) + PZPY
                    SMAT(INDI+4,INDJ+4) = SMAT(INDI+4,INDJ+4) + PZPZ

                 ENDIF



              ENDIF
           ENDIF
        ENDIF
     ENDDO

     !     IF (BASISI .EQ. "sp") INDI = INDI + 4
     !     IF (BASISI .EQ. "s") INDI = INDI + 1

  ENDDO

  !$OMP END PARALLEL DO

  DO I = 1, HDIM
     HDIAG(I) = H(I,I)
  ENDDO

  ! 
  ! ... and fill in the other half of the matrix
  !

  IF (BASISTYPE .EQ. "ORTHO") THEN

     DO I = 1, HDIM
        DO J = I, HDIM
           H(J,I) = H(I,J)
        ENDDO
     ENDDO

     !     NNZ = 0
     !     DO I = 1, HDIM
     !        DO J = 1, HDIM
     !           IF (ABS(H(J,I)) .GT. 1.0D-12) NNZ = NNZ + 1
     !        ENDDO
     !     ENDDO


     !     OPEN(UNIT=50, STATUS="UNKNOWN", FILE="myham.mtx")

     !     WRITE(50,'("%%MatrixMarket matrix coordinate real general")')
     !     WRITE(50,*) HDIM, HDIM, NNZ
     !     DO I = 1, HDIM
     !        DO J = 1, HDIM
     !           IF (ABS(H(I,J)) .GT. 1.0D-12)  WRITE(50,50) I, J, H(I,J)
     !        ENDDO
     !     ENDDO

     !50   FORMAT(2I8, G24.12)
     !     STOP

  ELSE

     DO I = 1, HDIM
        DO J = I, HDIM
           H(J,I) = H(I,J)
           SMAT(J,I) = SMAT(I,J)
        ENDDO
     ENDDO

     H0 = H

#ifdef PROGRESSON
     IF(LATTEINEXISTS)THEN
        CALL GENXBML
     ELSE
        CALL GENX
     ENDIF
#else
     CALL GENX
#endif

     IF (DEBUGON .EQ. 1) THEN

        OPEN(UNIT=30, STATUS="UNKNOWN", FILE="myS.dat")
        OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myH0.dat")

        PRINT*, "Caution - the Slater-Koster H and overlap matrices are being written to file"

        DO I = 1, HDIM
           WRITE(30,10) (SMAT(I,J), J = 1, HDIM)
           WRITE(31,10) (H0(I,J), J = 1, HDIM)
        ENDDO

        CLOSE(30)
        CLOSE(31)

     ENDIF

  ENDIF

  !  IF (BASISTYPE .EQ. "NONORTHO") CALL GENX

10 FORMAT(100F12.6)

  RETURN

END SUBROUTINE BLDNEWHS_SP
