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

SUBROUTINE FSPINNONO

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

  INTEGER :: I, J, K, L, M, N, KK, INDI, INDJ
  INTEGER :: LBRA, MBRA, LKET, MKET
  INTEGER :: PREVJ, NEWJ
  INTEGER :: PBCI, PBCJ, PBCK
  INTEGER :: BASISI(5), BASISJ(5), LBRAINC, LKETINC
  INTEGER :: SPININDI, SPININDJ
  REAL(LATTEPREC) :: ALPHA, BETA, PHI, RHO, COSBETA
  REAL(LATTEPREC) :: RIJ(3), DC(3)
  REAL(LATTEPREC) :: MAGR, MAGR2, MAGRP, MAGRP2, FTMP(3)
  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2
  REAL(LATTEPREC) :: MYDFDA, MYDFDB, MYDFDR, RCUTTB
  REAL(LATTEPREC) :: WSPINI, WSPINJ
  REAL(LATTEPREC), EXTERNAL :: DFDA, DFDB, DFDR
  LOGICAL PATH
  IF (EXISTERROR) RETURN

  FSSPIN = ZERO
  VIRSSPIN = ZERO

  !$OMP PARALLEL DO DEFAULT (NONE) &                                              
  !$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &                      
  !$OMP SHARED(CR, BOX, BO, RHOUP, RHODOWN, SPINON, NOINT, ATELE, ELE1, ELE2) &   
  !$OMP SHARED(BOND, OVERL, MATINDLIST, BASISTYPE) &                              
  !$OMP SHARED(DELTASPIN, WSS, WPP, WDD, WFF, SPININDLIST) &
  !$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &    
  !$OMP PRIVATE(RIJ, MAGR2, MAGR, MAGRP2, MAGRP, PATH, PHI, ALPHA, BETA, COSBETA, FTMP) &
  !$OMP PRIVATE(DC, LBRAINC, LBRA, MBRA, L, LKETINC, LKET, MKET, RHO) &           
  !$OMP PRIVATE(MYDFDA, MYDFDB, MYDFDR, RCUTTB, WSPINI, WSPINJ) &
  !$OMP PRIVATE(SPININDI, SPININDJ) &
  !$OMP REDUCTION(+:FSSPIN, VIRSSPIN)

  DO I = 1, NATS

     ! Build list of orbitals on atom I                                           
     SELECT CASE(BASIS(ELEMPOINTER(I)))

     CASE("s")
        BASISI(1) = 0
        BASISI(2) = -1
     CASE("p")
        BASISI(1) = 1
        BASISI(2) = -1
     CASE("d")
        BASISI(1) = 2
        BASISI(2) = -1
     CASE("f")
        BASISI(1) = 3
        BASISI(2) = -1
     CASE("sp")
        BASISI(1) = 0
        BASISI(2) = 1
        BASISI(3) = -1
     CASE("sd")
        BASISI(1) = 0
        BASISI(2) = 2
        BASISI(3) = -1
     CASE("sf")
        BASISI(1) = 0
        BASISI(2) = 3
        BASISI(3) = -1
     CASE("pd")
        BASISI(1) = 1
        BASISI(2) = 2
        BASISI(3) = -1
     CASE("pf")
        BASISI(1) = 1
        BASISI(2) = 3
        BASISI(3) = -1
     CASE("df")
        BASISI(1) = 2
        BASISI(2) = 3
        BASISI(3) = -1
     CASE("spd")
        BASISI(1) = 0
        BASISI(2) = 1
        BASISI(3) = 2
        BASISI(4) = -1
     CASE("spf")
        BASISI(1) = 0
        BASISI(2) = 1
        BASISI(3) = 3
        BASISI(4) = -1
     CASE("sdf")
        BASISI(1) = 0
        BASISI(2) = 2
        BASISI(3) = 3
        BASISI(4) = -1
     CASE("pdf")
        BASISI(1) = 1
        BASISI(2) = 2
        BASISI(3) = 3
        BASISI(4) = -1
     CASE("spdf")
        BASISI(1) = 0
        BASISI(2) = 1
        BASISI(3) = 2
        BASISI(4) = 3
        BASISI(5) = -1
     END SELECT

     INDI = MATINDLIST(I)
     SPININDI = SPININDLIST(I)

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

        RCUTTB = ZERO
        DO K = 1, NOINT

           IF ( (ATELE(I) .EQ. ELE1(K) .AND. &
                ATELE(J) .EQ. ELE2(K)) .OR. &
                (ATELE(J) .EQ. ELE1(K) .AND. &
                ATELE(I) .EQ. ELE2(K) )) THEN

              IF (BOND(8,K) .GT. RCUTTB ) RCUTTB = BOND(8,K)

              IF (BASISTYPE .EQ. "NONORTHO") THEN
                 IF (OVERL(8,K) .GT. RCUTTB ) RCUTTB = OVERL(8,K)
              ENDIF

           ENDIF

        ENDDO

        IF (MAGR2 .LT. RCUTTB*RCUTTB) THEN

           MAGR = SQRT(MAGR2)

           ! Build list of orbitals on atom J                                     

           SELECT CASE(BASIS(ELEMPOINTER(J)))
           CASE("s")
              BASISJ(1) = 0
              BASISJ(2) = -1
           CASE("p")
              BASISJ(1) = 1
              BASISJ(2) = -1
           CASE("d")
              BASISJ(1) = 2
              BASISJ(2) = -1
           CASE("f")
              BASISJ(1) = 3
              BASISJ(2) = -1
           CASE("sp")
              BASISJ(1) = 0
              BASISJ(2) = 1
              BASISJ(3) = -1
           CASE("sd")
              BASISJ(1) = 0
              BASISJ(2) = 2
              BASISJ(3) = -1
           CASE("sf")
              BASISJ(1) = 0
              BASISJ(2) = 3
              BASISJ(3) = -1
           CASE("pd")
              BASISJ(1) = 1
              BASISJ(2) = 2
              BASISJ(3) = -1
           CASE("pf")
              BASISJ(1) = 1
              BASISJ(2) = 3
              BASISJ(3) = -1
           CASE("df")
              BASISJ(1) = 2
              BASISJ(2) = 3
              BASISJ(3) = -1
           CASE("spd")
              BASISJ(1) = 0
              BASISJ(2) = 1
              BASISJ(3) = 2
              BASISJ(4) = -1
           CASE("spf")
              BASISJ(1) = 0
              BASISJ(2) = 1
              BASISJ(3) = 3
              BASISJ(4) = -1
           CASE("sdf")
              BASISJ(1) = 0
              BASISJ(2) = 2
              BASISJ(3) = 3
              BASISJ(4) = -1
           CASE("pdf")
              BASISJ(1) = 1
              BASISJ(2) = 2
              BASISJ(3) = 3
              BASISJ(4) = -1
           CASE("spdf")
              BASISJ(1) = 0
              BASISJ(2) = 1
              BASISJ(3) = 2
              BASISJ(4) = 3
              BASISJ(5) = -1
           END SELECT

           INDJ = MATINDLIST(J)
           SPININDJ = SPININDLIST(J)

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

           FTMP = ZERO
           K = INDI

           LBRAINC = 1
           DO WHILE (BASISI(LBRAINC) .NE. -1)

              LBRA = BASISI(LBRAINC)
              LBRAINC = LBRAINC + 1

              SELECT CASE(LBRA)
              CASE(0)
                 WSPINI = WSS(ELEMPOINTER(I))
              CASE(1)
                 WSPINI = WPP(ELEMPOINTER(I))
              CASE(2)
                 WSPINI = WDD(ELEMPOINTER(I))
              CASE(3)
                 WSPINI = WFF(ELEMPOINTER(I))
              END SELECT

              WSPINI = WSPINI*DELTASPIN(SPININDI + LBRAINC - 1)

              DO MBRA = -LBRA, LBRA

                 K = K + 1
                 L = INDJ

                 LKETINC = 1
                 DO WHILE (BASISJ(LKETINC) .NE. -1)

                    LKET = BASISJ(LKETINC)
                    LKETINC = LKETINC + 1

                    SELECT CASE(LKET)
                    CASE(0)
                       WSPINJ = WSS(ELEMPOINTER(J))
                    CASE(1)
                       WSPINJ = WPP(ELEMPOINTER(J))
                    CASE(2)
                       WSPINJ = WDD(ELEMPOINTER(J))
                    CASE(3)
                       WSPINJ = WFF(ELEMPOINTER(J))
                    END SELECT

                    WSPINJ = WSPINJ*DELTASPIN(SPININDJ + LKETINC - 1)

                    DO MKET = -LKET, LKET

                       L = L + 1

                       RHO = RHOUP(L, K) - RHODOWN(L, K)

                       IF (.NOT. PATH) THEN

                          ! Unroll loops and pre-compute                          

                          MYDFDA = DFDA(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ALPHA, COSBETA, "S")

                          MYDFDB = DFDB(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ALPHA, COSBETA, "S")

                          MYDFDR = DFDR(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ALPHA, COSBETA, "S")

                          MYDFDA = MYDFDA * (WSPINI + WSPINJ)
                          MYDFDB = MYDFDB * (WSPINI + WSPINJ)
                          MYDFDR = MYDFDR * (WSPINI + WSPINJ)

                          !                                                       
                          ! d/d_alpha                                             
                          !                                                       

                          FTMP(1) = FTMP(1) + RHO * &
                               (-RIJ(2) / MAGRP2 * MYDFDA)

                          FTMP(2) = FTMP(2) + RHO * &
                               (RIJ(1)/ MAGRP2 * MYDFDA)

                          !                                                       
                          ! d/d_beta                                              
                          !                                                       

                          FTMP(1) = FTMP(1) + RHO * &
                               (((((RIJ(3) * RIJ(1)) / &
                               MAGR2)) / MAGRP) * MYDFDB)

                          FTMP(2) = FTMP(2) + RHO * &
                               (((((RIJ(3) * RIJ(2)) / &
                               MAGR2)) / MAGRP) * MYDFDB)

                          FTMP(3) = FTMP(3) - RHO * &
                               (((ONE - ((RIJ(3) * RIJ(3)) / &
                               MAGR2)) / MAGRP) * MYDFDB)

                          !                                                       
                          ! d/dR                                                  
                          !                                                       

                          FTMP(1) = FTMP(1) - RHO * DC(1) * &
                               MYDFDR

                          FTMP(2) = FTMP(2) - RHO * DC(2) * &
                               MYDFDR

                          FTMP(3) = FTMP(3) - RHO * DC(3) * &
                               MYDFDR


                       ELSE

                          ! pathological configuration in which beta=0            
                          ! or pi => alpha undefined                              

                          ! fixed: MJC 12/17/13                                   

                          MYDFDB = DFDB(I, J, LBRA, LKET, &
                               MBRA, MKET, MAGR, ZERO, COSBETA, "S") / MAGR

                          MYDFDB = MYDFDB * (WSPINI + WSPINJ)

                          FTMP(1) = FTMP(1) - RHO * (COSBETA * MYDFDB)

                          MYDFDB = DFDB(I, J, LBRA, LKET, &
                               MBRA, MKET, MAGR, PI/TWO, COSBETA, "S") / MAGR

                          MYDFDB = MYDFDB * (WSPINI + WSPINJ)

                          FTMP(2) = FTMP(2) - RHO * (COSBETA * MYDFDB)

                          MYDFDR = DFDR(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ZERO, COSBETA, "S")

                          MYDFDR = MYDFDR * (WSPINI + WSPINJ)

                          FTMP(3) = FTMP(3) - RHO * COSBETA * MYDFDR

                       ENDIF

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           !           FTMP = FTMP * ( HUBBARDU(ELEMPOINTER(J))*DELTAQ(J) + COULOMBV(J) &
           !                +HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) + COULOMBV(I))


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


        ENDIF
     ENDDO

  ENDDO

  !$OMP END PARALLEL DO

  VIRSSPIN = VIRSSPIN/TWO

  !  PRINT*, FSSPIN(1,1)

  !  DO I = 1, NATS
  !     WRITE(6,10) I, FSSPIN(1,I), FSSPIN(2,I), FSSPIN(3,I)
  !  ENDDO

  !10 FORMAT(I4, 3F12.6)

  RETURN

END SUBROUTINE FSPINNONO
