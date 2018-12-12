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

SUBROUTINE SUMMARY

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE SPARSEARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE MYPRECISION
  USE DIAGARRAY
  USE VIRIALARRAY
  USE NONOARRAY

  IMPLICIT NONE

  INTEGER :: I, J, INDEX, NUMORB
  REAL(LATTEPREC) :: FDIRAC, SUMQ, ATOMSPIN
  IF (EXISTERROR) RETURN

  IF( VERBOSE < 0 ) RETURN

  CALL WRTCFGS(-999)

  OPEN(UNIT=24, STATUS="UNKNOWN", FILE="mylastLATTEcalc")

  IF (MDON .EQ. 1) THEN
     WRITE(6,'("# Molecular dynamics calculation")')
     WRITE(24,'("Molecular dynamics calculation")')
     IF (QITER .EQ. 0) THEN
        WRITE(6,'("# Using SCF-free Fast QMD")')
        WRITE(24,'("Using SCF-free Fast QMD")')
     ENDIF
  ENDIF


  IF (ELECTRO .EQ. 1) THEN
     WRITE(6,'("# Self-consistent charge transfer: on")')    
     WRITE(6,'("# SCF tolerance = ", G8.3)') ELEC_QTOL
     WRITE(24,'("Self-consistent charge transfer: on")')    
     WRITE(24,'("SCF tolerance = ", G8.3)') ELEC_QTOL
     IF (ELECMETH .EQ. 0) THEN
        WRITE(6,'("# Using Ewald summation")')
        WRITE(6,'("# Coulomb accuracy = ",G8.3)') COULACC
        WRITE(6,'("# Real space cut-off for Ewald sum = ", F6.2)') COULCUT
        WRITE(24,'("Using Ewald summation")')
        WRITE(24,'("Coulomb accuracy = ",G8.3)') COULACC
        WRITE(24,'("Real space cut-off for Ewald sum = ", F6.2)') COULCUT
     ELSEIF (ELECMETH .EQ. 1) THEN
        WRITE(6,'("# Real space electrostatics")')
        WRITE(6,'("# Cut-off tail applied between ", F6.2, F6.2)') &
             COULR1, COULCUT
        WRITE(24,'("Real space electrostatics")')
        WRITE(24,'("Cut-off tail applied between ", F6.2, F6.2)') &
             COULR1, COULCUT
     ENDIF
  ELSE
     WRITE(6,'("# Local charge neutrality: on")')    
     WRITE(6,'("# SCF tolerance = ", G8.3)') ELEC_QTOL
     WRITE(24,'("Local charge neutrality: on")')    
     WRITE(24,'("SCF tolerance = ", G8.3)') ELEC_QTOL
  ENDIF

  IF (SPINON .EQ. 1) THEN
     WRITE(6,'("# Spin-polarized calculation")')
     WRITE(24,'("Spin-polarized calculation")')
  ENDIF

  IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1) THEN
     WRITE(6,'("# Sparse matrix calculations for O(N)")')
     WRITE(24,'("Sparse matrix calculations for O(N)")')
     WRITE(6,'("# Numerical threshold = ",G8.3)') NUMTHRESH
     WRITE(24,'("Numerical threshold = ",G8.3)') NUMTHRESH
  ENDIF

  IF (CONTROL .EQ. 1) THEN
     WRITE(6,'("# Diagonalization")')
     WRITE(24,'("Diagonalization")')
  ELSEIF (CONTROL .EQ. 2) THEN
     WRITE(6,'("# SP2 purification at zero temperature")')
     WRITE(24,'("SP2 purification at zero temperature")')
  ELSEIF (CONTROL .EQ. 3) THEN
     WRITE(6,'("# Recursive expansion of the Fermi operator")')
     WRITE(24,'("Recursive expansion of the Fermi operator")')
  ELSEIF (CONTROL .EQ. 4) THEN
     WRITE(6,'("# SP2 purification at finite temperature")')
     WRITE(24,'("SP2 purification at finite temperature")')
  ELSEIF (CONTROL .EQ. 5) THEN
     WRITE(6,'("# SP2/Fermi method at finite temperature")')
     WRITE(24,'("SP2/Fermi method at finite temperature")')
  ENDIF

  !  CALL ALLOCATEDIAG

  IF (ALLOCATED(DIAG_IWORK) .EQV. .TRUE.) THEN
     WRITE(6,'("# Using xSYEVD for diaginalization - potentially fast and potentially unreliable")')
     WRITE(24,'("# Using xSYEVD for diaginalization - potentially fast and potentially unreliable")')
  ELSEIF (ALLOCATED(DIAG_IWORK) .EQV. .FALSE.) THEN
     WRITE(6,'("# Using xSYEV for diaginalization")')
     WRITE(24,'("# Using xSYEV for diaginalization")')
  ENDIF


  IF (CONTROL .NE. 2) THEN
     WRITE(6,'("# KBT (in eV) = ", F16.8)') KBT
     WRITE(24,'("KBT (in eV) = ", F16.8)') KBT
  ENDIF

  IF (CONTROL .EQ. 5) THEN
     WRITE(6,'("# Gershgorin: MAXEVAL, MINEVAL = ", 2F16.8)') MAXEVAL, MINEVAL
     WRITE(24,'("Gershgorin: MAXEVAL, MINEVAL = ", 2F16.8)') MAXEVAL, MINEVAL
  ENDIF

  IF (LATTEPREC .EQ. KIND(0.0D0)) THEN
     WRITE(6,'("# Double precision arithmetic")')
     WRITE(24,'("Double precision arithmetic")')
  ELSEIF (LATTEPREC .EQ. KIND(0.0)) THEN
     WRITE(6,'("# Single precision arithmetic")')
     WRITE(24,'("Single precision arithmetic")')
  ENDIF

  IF (ENTROPYKIND .EQ. 0) THEN
     WRITE(6,'("# Entropy set = 0")')
     WRITE(24,'("Entropy set = 0")')
  ELSEIF (ENTROPYKIND .EQ. 1) THEN
     WRITE(6,'("# Using exact ln form for entropy")')
     WRITE(24,'("Using exact ln form for entropy")')
  ELSEIF (ENTROPYKIND .EQ. 2) THEN
     WRITE(6,'("# Using the close-to-exact expansion of exact entropy (2)")')
     WRITE(24,'("Using the close-to-exact expansion of exact entropy (2)")')
  ELSEIF (ENTROPYKIND .EQ. 3) THEN
     WRITE(6,'("# Using 4th order approximation for entropy")')
     WRITE(24,'("Using 4th order approximation for entropy")')
  ELSEIF (ENTROPYKIND .EQ. 4) THEN
     WRITE(6,'("# Using 8th order approximation for entropy")')
     WRITE(24,'("Using 8th order approximation for entropy")') 
  ENDIF

  WRITE(6,'("# Tr[ rho*H ] = ", F16.8)') TRRHOH
  WRITE(6,'("# Pairwise energy = ", F16.8)') EREP
  WRITE(24,'("Tr[ rho*H ] = ", F16.8)') TRRHOH
  WRITE(24,'("Pairwise energy = ", F16.8)') EREP

  IF (ELECTRO .EQ. 1) THEN
     WRITE(6,'("# Coulombic + onsite E = ", F16.8)') ECOUL
     WRITE(24,'("Coulombic + onsite E = ", F16.8)') ECOUL
  ENDIF

  IF (CONTROL .NE. 2) THEN
     WRITE(6,'("# Electron entropy TS = ", F16.8)') ENTE
     WRITE(24,'("Electron entropy TS = ", F16.8)') ENTE
  ENDIF

  IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 .OR. CONTROL .EQ. 5) THEN
     WRITE(6,'("# Chemical potential = ", F16.8)') CHEMPOT
     WRITE(24,'("Chemical potential = ", F16.8)') CHEMPOT
  ENDIF

  IF (SPINON .EQ. 1) THEN
     WRITE(6,'("# Self-consistent spin energy = ", F16.8)') ESPIN
     WRITE(6,'("# Free atom spin energy = ", F16.8)') ESPIN_ZERO
     WRITE(24,'("Self-consistent spin energy = ", F16.8)') ESPIN
     WRITE(24,'("Free atom spin energy = ", F16.8)') ESPIN_ZERO
  ENDIF

  CALL GETPRESSURE

  WRITE(6,'("# Pressure (GPa) = ", F16.8)') PRESSURE
  WRITE(24,'("# Pressure (GPa) = ", F16.8)') PRESSURE

  IF (SPINON .EQ. 0) THEN
     WRITE(6,'("# Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL
     WRITE(6,'("")')
     WRITE(6,'("# FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE
     WRITE(6,'("")')
     WRITE(24,'("Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL
     WRITE(24,'("")')
     WRITE(24,'("FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE
     WRITE(24,'("")')
  ELSEIF (SPINON .EQ. 1) THEN
     WRITE(6,'("# Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL + &
          ESPIN 
     WRITE(6,'("")')
     WRITE(6,'("# FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE + &
          ESPIN 
     WRITE(6,'("")')
     WRITE(24,'("Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL + &
          ESPIN 
     WRITE(24,'("")')
     WRITE(24,'("FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE + &
          ESPIN 
     WRITE(24,'("")')
  ENDIF

  WRITE(6,60) "#checkP ", SYSVOL, TRRHOH + EREP - ECOUL - ENTE, -PRESSURE/TOGPA

60 FORMAT(A7, 3F18.6)

  ! Write the stress tensor

  WRITE(6,'("")')
  WRITE(6,'("# Stress tensor (GPa)")')
  WRITE(6,251) "# ", STRTEN(1), STRTEN(4), STRTEN(6)
  WRITE(6,251) "# ", STRTEN(4), STRTEN(2), STRTEN(5)
  WRITE(6,251) "# ", STRTEN(6), STRTEN(5), STRTEN(3)
  WRITE(6,'("")')
  WRITE(24,'("")')
  WRITE(24,'("# Stress tensor (GPa)")')
  WRITE(24,251) "# ",  STRTEN(1), STRTEN(4), STRTEN(6)
  WRITE(24,251) "# ", STRTEN(4), STRTEN(2), STRTEN(5)
  WRITE(24,251) "# ", STRTEN(6), STRTEN(5), STRTEN(3)
  WRITE(24,'("")')

250 FORMAT(3G20.9)
251 FORMAT(A2,1X,3G20.9)

  IF (KON .EQ. 0) THEN

     IF (SPINON .EQ. 0) THEN

        ! Analyse the eigenvalues only if we've used diagonalization
        ! to compute the density matrix

        IF (CONTROL .EQ. 1) THEN

           IF (BASISTYPE .EQ. "ORTHO") THEN

              WRITE(6,'("")')
              WRITE(6,'("# Eigenvalues")')
              WRITE(24,'("")')
              WRITE(24,'("Eigenvalues")')

              IF (KBT .GT. 0.00000001) THEN

                 DO I = 1, HDIM

                    FDIRAC = ONE/(ONE+EXP((EVALS(I) - CHEMPOT)/KBT))

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "*"
                       WRITE(24, 53) I, EVALS(I),  FDIRAC, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "o"
                       WRITE(24, 53) I, EVALS(I), FDIRAC, "o"

                    ENDIF

                 ENDDO

              ELSE

                 DO I = 1, HDIM

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), ONE, "*"
                       WRITE(24, 53) I, EVALS(I),  ONE, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), ZERO, "o"
                       WRITE(24, 53) I, EVALS(I), ZERO, "o"

                    ENDIF

                 ENDDO

              ENDIF

              WRITE(6,'("")')
              WRITE(6,'("# Band width = ", F12.6)') EVALS(HDIM) - EVALS(1)
              WRITE(6,'("")')
              WRITE(24,'("")')
              WRITE(24,'("Band width = ", F12.6)') EVALS(HDIM) - EVALS(1)
              WRITE(24,'("")')

              WRITE(6,'("")')
              WRITE(6,'("# HOMO-LUMO = ", F12.6)') EGAP
              WRITE(6,'("")')
              WRITE(24,'("")')
              WRITE(24,'("HOMO-LUMO = ", F12.6)') EGAP
              WRITE(24,'("")')


           ELSE

              ! Non-orthogonal basis - look at eigenvalus of H and XHX

              WRITE(6,'("")')
              WRITE(6,'("# Eigenvalues XHX")')
              WRITE(24,'("")')
              WRITE(24,'("Eigenvalues XHX")')

              IF (KBT .GT. 0.00000001) THEN

                 DO I = 1, HDIM

                    FDIRAC = ONE/(ONE+EXP((EVALS(I) - CHEMPOT)/KBT))

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "*"
                       WRITE(24, 53) I, EVALS(I),  FDIRAC, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "o"
                       WRITE(24, 53) I, EVALS(I), FDIRAC, "o"

                    ENDIF

                 ENDDO

                 !              WRITE(6,'("Eigenvectors of XHX")')
                 !              WRITE(24,'("Eigenvectors of XHX")')
                 ! Eigenvectors too 

                 !              DO I = 1, HDIM

                 !                 WRITE(6,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)
                 !                 WRITE(24,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)

                 !              ENDDO

              ELSE

                 DO I = 1, HDIM

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), ONE, "*"
                       WRITE(24, 53) I, EVALS(I),  ONE, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), ZERO, "o"
                       WRITE(24, 53) I, EVALS(I), ZERO, "o"

                    ENDIF

                 ENDDO

                 !              WRITE(6,'("Eigenvectors of XHX")')
                 !              WRITE(24,'("Eigenvectors of XHX")')
                 ! Eigenvectors too                                                                                 

                 !              DO I = 1, HDIM

                 !                 WRITE(6,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)
                 !                 WRITE(24,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)

                 !              ENDDO

              ENDIF

              ORTHOH = H

              CALL DIAGMYH

              WRITE(6,'("")')
              WRITE(6,'("# Eigenvalues H")')
              WRITE(24,'("")')
              WRITE(24,'("Eigenvalues H")')

              IF (KBT .GT. 0.00000001) THEN

                 DO I = 1, HDIM

                    FDIRAC = ONE/(ONE+EXP((EVALS(I) - CHEMPOT)/KBT))

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "*"
                       WRITE(24, 53) I, EVALS(I),  FDIRAC, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), FDIRAC, "o"
                       WRITE(24, 53) I, EVALS(I), FDIRAC, "o"

                    ENDIF

                 ENDDO

                 !              WRITE(6,'("Eigenvectors of H")')
                 !              WRITE(24,'("Eigenvectors of H")')
                 ! Eigenvectors too 

                 !             DO I = 1, HDIM

                 !                WRITE(6,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)
                 !                WRITE(24,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)

                 !             ENDDO


              ELSE

                 DO I = 1, HDIM

                    IF (EVALS(I) .LE. CHEMPOT) THEN

                       WRITE(6, 54) "# ", I, EVALS(I), ONE, "*"
                       WRITE(24, 53) I, EVALS(I),  ONE, "*"

                    ELSE 

                       WRITE(6, 54) "# ", I, EVALS(I), ZERO, "o"
                       WRITE(24, 53) I, EVALS(I), ZERO, "o"

                    ENDIF

                 ENDDO

                 !              WRITE(6,'("Eigenvectors of H")')
                 !              WRITE(24,'("Eigenvectors of H")')
                 ! Eigenvectors too 

                 !             DO I = 1, HDIM

                 !                WRITE(6,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)
                 !                WRITE(24,'(100F12.6)') (EVECS(I,J), J = 1, HDIM)

                 !             ENDDO


              ENDIF

           ENDIF

           !        CALL GENDIAG


           CALL DEALLOCATEDIAG

        ENDIF

53      FORMAT(I6, 2X, F14.8, 1X, G18.12, 1X,  A1)
54      FORMAT(A2, I6, 2X, F14.8, 1X, G18.12, 1X, A1)

        !     CALL DEALLOCATEDIAG

     ELSEIF (SPINON .EQ. 1 .AND. CONTROL .EQ. 1) THEN

        WRITE(6,'("")')
        WRITE(6,'("# Eigenvalues")')
        WRITE(6,'("#        Up            :           Down")')
        WRITE(24,'("")')
        WRITE(24,'("Eigenvalues")')
        WRITE(24,'("        Up            :           Down")')
        DO I = 1, HDIM
           WRITE(6, 152) "# ", I, UPEVALS(I), DOWNEVALS(I)
           WRITE(24, 52) I, UPEVALS(I), DOWNEVALS(I)
        ENDDO

52      FORMAT(I6, 2X, F14.8, 4X, F14.8)
152     FORMAT(A2, I6, 2X, F14.8, 4X, F14.8)

        CALL DEALLOCATEDIAG

     ENDIF

  ENDIF
  !  IF (ELECTRO .EQ. 1) THEN

  WRITE(6,'("# Partial charges")')
  WRITE(6,'("#  Atom    Type      Charge (e)")')
  WRITE(24,'("Partial charges")')
  WRITE(24,'("  Atom    Type      Charge (e)")')
  SUMQ = ZERO
  DO I = 1, NATS
     WRITE(6,150) "# ", I, ATELE(I), -DELTAQ(I)
     WRITE(24,50) I, ATELE(I), -DELTAQ(I)
     SUMQ = SUMQ + DELTAQ(I)
  ENDDO

  WRITE(6,'("# Mulliken occupancies")')
  WRITE(6,'("#  Atom    Type     Free atom     Self-consistent")')
  WRITE(24,'("Mulliken occupancies")')
  WRITE(24,'("  Atom    Type     Free atom     Self-consistent")')

  DO I = 1, NATS
     WRITE(6,155) "# ", I, ATELE(I), ATOCC(ELEMPOINTER(I)), MYCHARGE(I)
     WRITE(24,55) I, ATELE(I), ATOCC(ELEMPOINTER(I)), MYCHARGE(I)
  ENDDO


  WRITE(6,'("# Sum of partial charges =", G16.8)') SUMQ
  WRITE(24,'(" Sum of partial charges =", G16.8)') SUMQ

50 FORMAT(I6, 4X, A2, 2X, F11.8)
150 FORMAT(A2, I6, 4X, A2, 2X, F11.8)

55 FORMAT(I6, 4X, A2, 7X, F12.6, 5X, F12.6)
155 FORMAT(A2, I6, 4X, A2, 7X, F12.6, 5X, F12.6)

  !  ENDIF

  IF (SPINON .EQ. 1) THEN

     WRITE(6,'("")')
     WRITE(6,'("# Orbital spin densities")')
     WRITE(6,'("# Orbital index         Spin density")')
     WRITE(24,'("")')
     WRITE(24,'("Orbital spin densities")')
     WRITE(24,'("Orbital index         Spin density")')
     DO I = 1, DELTADIM
        WRITE(6,151) "# ", I, DELTASPIN(I)
        WRITE(24,51) I, DELTASPIN(I)
     ENDDO

     WRITE(6,'("")')
     WRITE(6,'("# Per atom spin densities")')
     WRITE(6,'("#   Atom   Type       Spin density")')
     WRITE(24,'("")')
     WRITE(24,'(" Per atom spin densities")')
     WRITE(24,'("   Atom   Type       Spin density")')

     INDEX = 0
     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           NUMORB = 1

        CASE("p")

           NUMORB = 1

        CASE("d")

           NUMORB = 1

        CASE("f")

           NUMORB = 1

        CASE("sp")

           NUMORB = 2

        CASE("sd")

           NUMORB = 2

        CASE("sf")

           NUMORB = 2

        CASE("pd")

           NUMORB = 2

        CASE("pf")

           NUMORB = 2

        CASE("df")

           NUMORB = 2

        CASE("spd")

           NUMORB = 3

        CASE("spf")

           NUMORB = 3

        CASE("sdf")

           NUMORB = 3

        CASE("pdf")

           NUMORB = 3

        CASE("spdf") 

           NUMORB = 4

        END SELECT

        ATOMSPIN = ZERO
        DO J = 1, NUMORB
           INDEX = INDEX + 1
           ATOMSPIN = ATOMSPIN + DELTASPIN(INDEX)
        ENDDO


        WRITE(6,156) "# ", I, ATELE(I), ATOMSPIN
        WRITE(24,56) I, ATELE(I), ATOMSPIN

     ENDDO


156  FORMAT(A2,1X,I6,3X,A2,6X,F12.6)
56   FORMAT(I6,1X,A2,1X,F12.6) 



51   FORMAT(I6, 16X,F14.8)
151  FORMAT(A2, I6, 16X,F14.8)

  ENDIF

  CLOSE(24)

  !  CALL WRTCFGS(-999)


  RETURN

END SUBROUTINE SUMMARY
