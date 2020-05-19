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
  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2, MYANGFACTOR
  REAL(LATTEPREC) :: AMMBRA, WIGLBRAMBRA
  REAL(LATTEPREC) :: MYBONDINT(0:3), MYOVERLAPINT(0:3)
  REAL(LATTEPREC), EXTERNAL :: ANGFACTOR, UNIVSCALE
  IF (EXISTERROR) RETURN


  H = ZERO

  IF (BASISTYPE .EQ. "NONORTHO") THEN

     SMAT = ZERO
     DO I = 1, HDIM
        SMAT(I,I) = ONE
     ENDDO

  ENDIF

  INDEX = 0

  ! Build diagonal elements (pre-calculated)

  DO I = 1, HDIM
     
     H(I, I) = H_ONSITE(I)

  ENDDO


!$OMP PARALLEL DO DEFAULT (NONE) &
!$OMP SHARED(NATS, BASIS, ELEMPOINTER, TOTNEBTB, NEBTB) &
!$OMP SHARED(CR, BOX, H, SMAT) &
!$OMP SHARED(HCUT, SCUT, MATINDLIST, BASISTYPE, ORBITAL_LIST, CUTOFF_LIST) &
!$OMP PRIVATE(I, J, K, NEWJ, BASISI, BASISJ, INDI, INDJ, PBCI, PBCJ, PBCK) &
!$OMP PRIVATE(RIJ, MAGR2, MAGR, MAGRP, PHI, ALPHA, BETA, COSBETA) &
!$OMP PRIVATE(LBRAINC, LBRA, MBRA, L, LKETINC, LKET, MKET) &
!$OMP PRIVATE(RCUTTB, IBRA, IKET, MYANGFACTOR, MP) &
!$OMP PRIVATE(MYBONDINT, MYOVERLAPINT)

  ! open loop over atoms I in system
  DO I = 1, NATS

     ! Build the lists of orbitals on each atom

     BASISI(:) = ORBITAL_LIST(:,I)

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
        
        RCUTTB = CUTOFF_LIST(J,I)

        IF (MAGR2 .LT. RCUTTB*RCUTTB) THEN
           
           MAGR = SQRT(MAGR2)

           BASISJ(:) = ORBITAL_LIST(:,J)
           
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
                 
                 IKET = INDJ + 1
                 
                 LKETINC = 1
                 DO WHILE (BASISJ(LKETINC) .NE. -1)
                    
                    LKET = BASISJ(LKETINC)
                    LKETINC = LKETINC + 1
                    
                    ! Precompute the integrals outide the MKET loop

                    DO MP = 0, MIN(LBRA, LKET)
                       
                       MYBONDINT(MP) = UNIVSCALE(I, J, LBRA, LKET, MP, MAGR, "H")
                       
                       IF (BASISTYPE .EQ. "NONORTHO") &
                            MYOVERLAPINT(MP) = UNIVSCALE(I, J, LBRA, LKET, MP, MAGR, "S")

                    ENDDO

                    
                    DO MKET = -LKET, LKET
                       
                       ! This is the sigma bonds (mp = 0)
                       
                       ! Hamiltonian build
                       
                       ! Pre-compute the angular part so we can use it
                          ! again later if we're building the S matrix too
                       
                       MYANGFACTOR = ANGFACTOR(LBRA, LKET, MBRA, MKET, 0, ALPHA, COSBETA)

                       H(IBRA, IKET) = H(IBRA, IKET) + MYANGFACTOR * &
                            MYBONDINT(0)
                       
                       ! Overlap matrix build
                       
                       IF (BASISTYPE .EQ. "NONORTHO") THEN
                          
                          SMAT(IBRA, IKET) = SMAT(IBRA, IKET) + MYANGFACTOR * &
                               MYOVERLAPINT(0)
                          
                       ENDIF
                       
                       ! everything else
                       
                       DO MP = 1, MIN(LBRA, LKET)

                          MYANGFACTOR = ANGFACTOR(LBRA, LKET, MBRA, MKET, MP, ALPHA, COSBETA)
                          
                          H(IBRA, IKET) = H(IBRA, IKET) + MYANGFACTOR * &
                               MYBONDINT(MP)
                          
                          IF (BASISTYPE .EQ. "NONORTHO") THEN
                             
                             SMAT(IBRA, IKET) = SMAT(IBRA, IKET) + &
                                  MYANGFACTOR * MYOVERLAPINT(MP)
                             
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
