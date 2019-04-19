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

SUBROUTINE BLDNEWHS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NEBLISTARRAY
  USE XBOARRAY
  USE NONOARRAY
  USE UNIVARRAY
  USE MYPRECISION

#ifdef PROGRESSON
  USE GENXPROGRESS
#endif

  IMPLICIT NONE

  INTEGER :: I, J, NEWJ, K, L, II, JJ, KK, MM, MP, NN, SUBI
  INTEGER :: IBRA, IKET, LBRA, LKET, MBRA, MKET
  INTEGER :: INDEX, INDI, INDJ
  INTEGER :: SWITCH, PREVJ
  INTEGER :: PBCI, PBCJ, PBCK
  INTEGER :: BASISI(5), BASISJ(5), LBRAINC, LKETINC
  REAL(LATTEPREC) :: ALPHA, BETA, COSBETA, PHI, TMP, PERM
  REAL(LATTEPREC) :: RIJ(3), MAGR2, MAGR, MAGRP, RCUTTB
  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2
  REAL(LATTEPREC) :: ANGFACTOR, AMMBRA, WIGLBRAMBRA
  !  REAL(LATTEPREC), ALLOCATABLE :: SAVEHDIAG(:)
  REAL(LATTEPREC), EXTERNAL :: UNIVSCALE, WIGNERD, SLMMP, TLMMP, AM, BM
  IF (EXISTERROR) RETURN


  H = ZERO

  IF (BASISTYPE .EQ. "NONORTHO") THEN

     SMAT = ZERO
     DO I = 1, HDIM
        SMAT(I,I) = ONE
     ENDDO

  ENDIF

  INDEX = 0

  ! Build diagonal elements
  DO I = 1, NATS

     K = ELEMPOINTER(I)

     SELECT CASE(BASIS(K))

     CASE("s")

        INDEX = INDEX + 1
        H(INDEX, INDEX) = HES(K)

     CASE("p")

        DO SUBI = 1, 3
           INDEX = INDEX + 1
           H(INDEX,INDEX) = HEP(K)
        ENDDO

     CASE("d")

        DO SUBI = 1, 5
           INDEX = INDEX + 1
           H(INDEX,INDEX) = HED(K)
        ENDDO

     CASE("f")

        DO SUBI = 1, 7
           INDEX = INDEX + 1
           H(INDEX,INDEX) = HEF(K)
        ENDDO

     CASE("sp")

        DO SUBI = 1, 4

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSE
              H(INDEX,INDEX) = HEP(K)
           ENDIF

        ENDDO

     CASE("sd")

        DO SUBI = 1, 6

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSE
              H(INDEX,INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("sf")

        DO SUBI = 1, 8

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("pd")

        DO SUBI = 1, 8

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H(INDEX,INDEX) = HEP(K)
           ELSE
              H(INDEX,INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("pf")

        DO SUBI = 1, 10

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H(INDEX,INDEX) = HEP(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("df")

        DO SUBI = 1, 12

           INDEX = INDEX + 1
           IF (SUBI .LE. 5) THEN
              H(INDEX,INDEX) = HED(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("spd")

        DO SUBI = 1, 9

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H(INDEX,INDEX) = HEP(K)
           ELSE
              H(INDEX,INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("spf")

        DO SUBI = 1, 11

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H(INDEX,INDEX) = HEP(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("sdf")

        DO SUBI = 1, 13

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX,INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 6) THEN
              H(INDEX,INDEX) = HED(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO


     CASE("pdf")

        DO SUBI = 1, 15

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H(INDEX,INDEX) = HEP(K)
           ELSEIF (SUBI .GT. 3 .AND. SUBI .LE. 8) THEN
              H(INDEX,INDEX) = HED(K)
           ELSE
              H(INDEX,INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("spdf")

        DO SUBI = 1, 16

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H(INDEX, INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H(INDEX, INDEX) = HEP(K)
           ELSEIF (SUBI .GT. 4 .AND. SUBI .LE. 9) THEN
              H(INDEX, INDEX) = HED(K)
           ELSE
              H(INDEX, INDEX) = HEF(K)
           ENDIF

        ENDDO

     END SELECT

  ENDDO

!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &
!$OMP SHARED(CR, BOX, H, SMAT, NOINT, ATELE, ELE1, ELE2) &
!$OMP SHARED(HCUT, SCUT, MATINDLIST, BASISTYPE) &
!$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &
!$OMP PRIVATE(RIJ, MAGR2, MAGR, MAGRP, PHI, ALPHA, BETA, COSBETA) &
!$OMP PRIVATE(LBRAINC, LBRA, MBRA, L, LKETINC, LKET, MKET) &
!$OMP PRIVATE(RCUTTB, IBRA, IKET, AMMBRA, WIGLBRAMBRA, ANGFACTOR, MP)

  ! open loop over atoms I in system
  DO I = 1, NATS

     ! Build the lists of orbitals on each atom

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

     ! open loop over neighbors J of atom I
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
              
              IF (HCUT(K) .GT. RCUTTB ) RCUTTB = HCUT(K)
              
              IF (BASISTYPE .EQ. "NONORTHO") THEN
                 IF (SCUT(K) .GT. RCUTTB ) RCUTTB = SCUT(K)
              ENDIF
              
           ENDIF
           
        ENDDO
        
        IF (MAGR2 .LT. RCUTTB*RCUTTB) THEN
           
           MAGR = SQRT(MAGR2)
           
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
           
           MAGRP = SQRT(RIJ(1) * RIJ(1) + RIJ(2) * RIJ(2))
           
           ! transform to system in which z-axis is aligned with RIJ,
           IF (ABS(RIJ(1)) .GT. 1.0E-12) THEN
              
              IF (RIJ(1) .GT. ZERO .AND. RIJ(2) .GE. ZERO) THEN
                 PHI = ZERO
              ELSEIF (RIJ(1) .GT. ZERO .AND. RIJ(2) .LT. ZERO) THEN
                 PHI = TWO * PI
              ELSE
                 PHI = PI
              ENDIF
              ALPHA = ATAN(RIJ(2) / RIJ(1)) + PHI
              
           ELSEIF (ABS(RIJ(2)) .GT. 1.0E-12) THEN
              
              IF (RIJ(2) .GT. 1.0E-12) THEN
                 ALPHA = PI / TWO
              ELSE
                 ALPHA = THREE * PI / TWO
              ENDIF
              
           ELSE
              ! pathological case: beta=0 and alpha undefined, but
              ! this doesn't matter for matrix elements
              
              ALPHA = ZERO
              
           ENDIF
           
           COSBETA = RIJ(3)/MAGR
           BETA = ACOS(RIJ(3) / MAGR)
           
           ! Build matrix elements using eqns (1)-(9) in PRB 72 165107
           
           ! The loops over LBRA and LKET need to take into account
           ! the orbitals assigned to each atom, e.g., sd rather than
           ! spd...
           
           IBRA = INDI + 1
           
           LBRAINC = 1
           DO WHILE (BASISI(LBRAINC) .NE. -1)
              
              LBRA = BASISI(LBRAINC)
              LBRAINC = LBRAINC + 1
              
              DO MBRA = -LBRA, LBRA
                 
                 ! We can calculate these two outside the
                 ! MKET loop...
                 
                 AMMBRA = AM(MBRA, ALPHA)
                 WIGLBRAMBRA = WIGNERD(LBRA, ABS(MBRA), 0, COSBETA)
                 
                 IKET = INDJ + 1
                 
                 LKETINC = 1
                 DO WHILE (BASISJ(LKETINC) .NE. -1)
                    
                    LKET = BASISJ(LKETINC)
                    LKETINC = LKETINC + 1
                    
                    DO MKET = -LKET, LKET
                       
                       ! This is the sigma bonds (mp = 0)
                       
                       ! Hamiltonian build
                       
                       ! Pre-compute the angular part so we can use it
                          ! again later if we're building the S matrix too
                       
                       ANGFACTOR = TWO * AMMBRA * &
                            AM(MKET, ALPHA) * &
                            WIGLBRAMBRA  * &
                            WIGNERD(LKET, ABS(MKET), 0, COSBETA)
                       
                       H(IBRA, IKET) = H(IBRA, IKET) + ANGFACTOR * &
                            UNIVSCALE(I, J, LBRA, LKET, 0, MAGR, "H")
                       
                       ! Overlap matrix build
                       
                       IF (BASISTYPE .EQ. "NONORTHO") THEN
                          
                          SMAT(IBRA, IKET) = SMAT(IBRA, IKET) + ANGFACTOR * &
                               UNIVSCALE(I, J, LBRA, LKET, 0, MAGR, "S")
                          
                       ENDIF
                       
                       ! everything else
                       
                       DO MP = 1, MIN(LBRA, LKET)
                          
                          ANGFACTOR = SLMMP(LBRA, MBRA, MP, ALPHA, COSBETA) * &
                               SLMMP(LKET, MKET, MP, ALPHA, COSBETA) + &
                               TLMMP(LBRA, MBRA, MP, ALPHA, COSBETA) * &
                               TLMMP(LKET, MKET, MP, ALPHA, COSBETA)
                          
                          H(IBRA, IKET) = H(IBRA, IKET) + ANGFACTOR * &
                               UNIVSCALE(I, J, LBRA, LKET, MP, MAGR, "H")
                          
                          IF (BASISTYPE .EQ. "NONORTHO") THEN
                             
                             SMAT(IBRA, IKET) = SMAT(IBRA, IKET) + &
                                  ANGFACTOR * &
                                  UNIVSCALE(I, J, LBRA, LKET, MP, MAGR, "S")
                             
                          ENDIF

                       ENDDO
                       
                       
                       IKET = IKET + 1

                    ENDDO
                    
                 ENDDO
                 
                 IBRA = IBRA + 1
                 
              ENDDO
           ENDDO
        ENDIF
     ENDDO

  ENDDO

!$OMP END PARALLEL DO

  DO I = 1, HDIM
     HDIAG(I) = H(I,I)
  ENDDO


  IF (BASISTYPE .EQ. "NONORTHO") THEN

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

10      FORMAT(100F12.6)

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE BLDNEWHS
