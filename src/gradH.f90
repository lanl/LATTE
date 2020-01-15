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

SUBROUTINE GRADH

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NEBLISTARRAY
  USE UNIVARRAY
  USE SPINARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, L, M, N, KK, INDI, INDJ
  INTEGER :: LBRA, MBRA, LKET, MKET
  INTEGER :: PREVJ, NEWJ, MP
  INTEGER :: PBCI, PBCJ, PBCK
  INTEGER :: BASISI(5), BASISJ(5), LBRAINC, LKETINC
  REAL(LATTEPREC) :: ALPHA, BETA, PHI, RHO, COSBETA
  REAL(LATTEPREC) :: RIJ(3), DC(3)
  REAL(LATTEPREC) :: MAGR, MAGR2, MAGRP, MAGRP2, FTMP(3)
  REAL(LATTEPREC) :: MYSLMMPL2(10), MYSLMMPL1(10), MYTMMPL2(10), MYTMMPL1(10), MYUNIVSCALE(10) 
  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2
  REAL(LATTEPREC) :: MYDFDA, MYDFDB, MYDFDR, RCUTTB
  REAL(LATTEPREC) :: MYDFDA1, MYDFDB1, MYDFDR1
  REAL(LATTEPREC), EXTERNAL :: DFDAPREC, DFDBPREC, DFDRPREC
  REAL(LATTEPREC), EXTERNAL :: DFDA, DFDB, DFDR
  REAL(LATTEPREC), EXTERNAL :: SLMMP, TLMMP, UNIVSCALE
  LOGICAL PATH
  IF (EXISTERROR) RETURN

  F = ZERO
  VIRBOND = ZERO

!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &
!$OMP SHARED(CR, BOX, BO, RHOUP, RHODOWN, SPINON) &
!$OMP SHARED(HCUT, SCUT, MATINDLIST, BASISTYPE, ORBITAL_LIST, CUTOFF_LIST) &
!$OMP SHARED(F) &
!$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &
!$OMP PRIVATE(RIJ, MAGR2, MAGR, MAGRP2, MAGRP, PATH, PHI, ALPHA, BETA, COSBETA, FTMP) &
!$OMP PRIVATE(DC, LBRAINC, LBRA, MBRA, L, LKETINC, LKET, MKET, RHO) &
!$OMP PRIVATE(MYDFDA, MYDFDB, MYDFDR, RCUTTB) &
!$OMP REDUCTION(+:VIRBOND)        

  DO I = 1, NATS

     ! Build list of orbitals on atom I

     BASISI(:) = ORBITAL_LIST(:,I)

     INDI = MATINDLIST(I)

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

           !        PRINT*, ALPHA, BETA

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
                       END SELECT

                       IF (.NOT. PATH) THEN

                          ! Unroll loops and pre-compute

                          CALL DFDX(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ALPHA, COSBETA, "H", &
                               MYDFDA, MYDFDB, MYDFDR)

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

                          CALL DFDX(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, ZERO, COSBETA, "H", &
                               MYDFDA, MYDFDB, MYDFDR)

                          MYDFDB = MYDFDB/MAGR

                          FTMP(1) = FTMP(1) - RHO * (COSBETA * MYDFDB)

                          FTMP(3) = FTMP(3) - RHO * COSBETA * MYDFDR 

                          
                          CALL DFDX(I, J, LBRA, LKET, MBRA, &
                               MKET, MAGR, PI/TWO, COSBETA, "H", &
                               MYDFDA, MYDFDB, MYDFDR)

                          MYDFDB = MYDFDB/MAGR

                          FTMP(2) = FTMP(2) - RHO * (COSBETA * MYDFDB)


                       ENDIF

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO

           F(1,I) = F(1,I) + FTMP(1)
           F(2,I) = F(2,I) + FTMP(2)
           F(3,I) = F(3,I) + FTMP(3)

           VIRBOND(1) = VIRBOND(1) + RIJ(1) * FTMP(1)
           VIRBOND(2) = VIRBOND(2) + RIJ(2) * FTMP(2)
           VIRBOND(3) = VIRBOND(3) + RIJ(3) * FTMP(3)
           VIRBOND(4) = VIRBOND(4) + RIJ(1) * FTMP(2)
           VIRBOND(5) = VIRBOND(5) + RIJ(2) * FTMP(3)
           VIRBOND(6) = VIRBOND(6) + RIJ(3) * FTMP(1)

        ENDIF

     ENDDO

     !     INDI = INDI + NORBI

  ENDDO

!$OMP END PARALLEL DO

  RETURN

END SUBROUTINE GRADH
