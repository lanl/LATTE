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

SUBROUTINE PANIC

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE SPARSEARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE MYPRECISION
  USE DIAGARRAY

  IMPLICIT NONE

  INTEGER :: I
  IF (EXISTERROR) RETURN


  IF (MDON .EQ. 1) THEN

     !
     ! Write a .cfg file with our problematic configuration
     !

     CALL WRTCFGS(-1)

  ELSE

     CALL WRTCFGS(-999)

  ENDIF

  OPEN(UNIT=24, STATUS="UNKNOWN", FILE="myPANICfile.dat")

  WRITE(6,'("# PANIC: SOMETHING BAD HAS HAPPENED! Check myPANICfile.dat")')
  WRITE(24,'("PANIC: SOMETHING BAD HAS HAPPENED! SOME CLUES TO DIAGNOSE THE PROBLEM FOLLOW")')

  IF (MDON .EQ. 1) THEN
     WRITE(24,'("Molecular dynamics calculation")')
  ENDIF


  IF (ELECTRO .EQ. 1) THEN
     WRITE(24,'("Self-consistent charge transfer: on")')
     WRITE(24,'("SCF tolerance = ", G8.3)') ELEC_QTOL
     IF (ELECMETH .EQ. 0) THEN
        WRITE(24,'("Using Ewald summation")')
        WRITE(24,'("Coulomb accuracy = ",G8.3)') COULACC
        WRITE(24,'("Real space cut-off for Ewald sum = ", F6.2)') COULCUT
     ELSEIF (ELECMETH .EQ. 1) THEN
        WRITE(24,'("Real space electrostatics")')
        WRITE(24,'("Cut-off tail applied between ", F6.2, F6.2)') &
             COULR1, COULCUT
     ENDIF
  ENDIF

  IF (SPINON .EQ. 1) THEN
     WRITE(24,'("Spin-polarized calculation")')
  ENDIF

  IF (SPARSEON .EQ. 1) THEN
     WRITE(24,'("Sparse matrix calculation")')
     WRITE(24,'("Numerical threshold = ",G8.3)') NUMTHRESH
  ENDIF



  IF (CONTROL .EQ. 1) THEN
     WRITE(24,'("Diagonalization")')
  ELSEIF (CONTROL .EQ. 2) THEN
     WRITE(24,'("SP2 purification at zero temperature")')
  ELSEIF (CONTROL .EQ. 3) THEN
     WRITE(24,'("Recursive expansion of the Fermi operator")')
  ELSEIF (CONTROL .EQ. 4) THEN
     WRITE(24,'("SP2 purification at finite temperature")')
  ELSEIF (CONTROL .EQ. 5) THEN
     WRITE(24,'("SP2/Fermi method at finite temperature")')
  ENDIF

  IF (CONTROL .NE. 2) THEN
     WRITE(24,'("KBT (in eV) = ", F16.8)') KBT
  ENDIF

  IF (CONTROL .EQ. 5) THEN
     WRITE(24,'("Gershgorin: MAXEVAL, MINEVAL = ", 2F16.8)') MAXEVAL, MINEVAL
  ENDIF

  IF (LATTEPREC .EQ. KIND(0.0D0)) THEN
     WRITE(24,'("Double precision arithmetic")')
  ELSEIF (LATTEPREC .EQ. KIND(0.0)) THEN
     WRITE(24,'("Single precision arithmetic")')
  ENDIF

  IF (ENTROPYKIND .EQ. 0) THEN
     WRITE(24,'("Entropy set = 0")')
  ELSEIF (ENTROPYKIND .EQ. 1) THEN
     WRITE(24,'("Using exact ln form for entropy")')
  ELSEIF (ENTROPYKIND .EQ. 2) THEN
     WRITE(24,'("Using the close-to-exact expansion of exact entropy (2)")')
  ELSEIF (ENTROPYKIND .EQ. 3) THEN
     WRITE(24,'("Using 4th order approximation for entropy")')
  ELSEIF (ENTROPYKIND .EQ. 4) THEN
     WRITE(24,'("Using 8th order approximation for entropy")')
  ENDIF

  WRITE(24,'("Tr[ rho*H ] = ", F16.8)') TRRHOH
  WRITE(24,'("Pairwise energy = ", F16.8)') EREP

  IF (ELECTRO .EQ. 1) THEN
     WRITE(24,'("Coulombic + onsite E = ", F16.8)') ECOUL
  ENDIF

  IF (CONTROL .NE. 2) THEN
     WRITE(24,'("Electron entropy TS = ", F16.8)') ENTE
  ENDIF

  IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 .OR. CONTROL .EQ. 5) THEN
     WRITE(24,'("Chemical potential = ", F16.8)') CHEMPOT
  ENDIF

  IF (SPINON .EQ. 1) THEN
     WRITE(24,'("Self-consistent spin energy = ", F16.8)') ESPIN
     WRITE(24,'("Free atom spin energy = ", F16.8)') ESPIN_ZERO
  ENDIF

  IF (SPINON .EQ. 0) THEN
     WRITE(24,'("Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL
     WRITE(24,'("")')
     WRITE(24,'("FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE
     WRITE(24,'("")')
  ELSEIF (SPINON .EQ. 1) THEN
     WRITE(24,'("Total energy (zero K) = ", F16.8)') TRRHOH + EREP - ECOUL + &
          ESPIN - ESPIN_ZERO
     WRITE(24,'("")')
     WRITE(24,'("FREE ENERGY = ", F16.8)') TRRHOH + EREP - ECOUL - ENTE + &
          ESPIN - ESPIN_ZERO
     WRITE(24,'("")')
  ENDIF

  IF (ELECTRO .EQ. 1) THEN

     WRITE(24,'("Partial charges")')
     WRITE(24,'("  Atom         Charge")')
     DO I = 1, NATS
        WRITE(24,50) I, DELTAQ(I)
     ENDDO

50   FORMAT(I6, 9X, F11.8)

  ENDIF

  IF (SPINON .EQ. 1) THEN

     WRITE(24,'("")')
     WRITE(24,'("Orbital spin densitites")')
     WRITE(24,'("Orbital index         Spin density")')
     DO I = 1, DELTADIM
        WRITE(24,51) I, DELTASPIN(I)
     ENDDO

51   FORMAT(I6, 16X,F14.8)

     IF (CONTROL .NE. 1) THEN
        CALL ALLOCATEDIAG
     ENDIF

     CALL DIAGMYH

     WRITE(24,'("")')
     WRITE(24,'("Eigenvalues")')
     WRITE(24,'("        Up            :           Down")')
     DO I = 1, HDIM
        WRITE(24, 52) I, UPEVALS(I), DOWNEVALS(I)
     ENDDO

52   FORMAT(I6, 2X, F14.8, 4X, F14.8)

     CALL DEALLOCATEDIAG

  ENDIF


  WRITE(24,'("")')
  WRITE(24,'("Coordinates (in .xyz format!)")')
  DO I = 1, NATS
     WRITE(24,54) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
  ENDDO

53 FORMAT(I6, 1X, 3F18.9, 1X, A2)
54 FORMAT(A2, 1X, 3G18.9)

  IF (MDON .EQ. 1) THEN

     WRITE(24,'("")')
     WRITE(24,'("Velocities")')
     DO I = 1, NATS
        WRITE(24,53) I, V(1,I), V(2,I), V(3,I), ATELE(I)
     ENDDO

  ENDIF

  WRITE(24,'("")')
  WRITE(24,'("Forces")')
  DO I = 1, NATS
     WRITE(24,53) I, FTOT(1,I), FTOT(2,I), FTOT(3,I), ATELE(I)
  ENDDO

  WRITE(24,'("")')
  WRITE(24,'("Band structure force")')
  DO I = 1, NATS
     WRITE(24,53) I, TWO*F(1,I), TWO*F(2,I), TWO*F(3,I), ATELE(I)
  ENDDO

  WRITE(24,'("")')
  WRITE(24,'("Pair potential force")')
  DO I = 1, NATS
     WRITE(24,53) I, FPP(1,I), FPP(2,I), FPP(3,I), ATELE(I)
  ENDDO

  IF (ELECTRO .EQ. 1) THEN

     WRITE(24,'("")')
     WRITE(24,'("Coulomb force")')
     DO I = 1, NATS
        WRITE(24,53) I, FCOUL(1,I), FCOUL(2,I), FCOUL(3,I), ATELE(I)
     ENDDO

  ENDIF

  CLOSE(24)

  RETURN

END SUBROUTINE PANIC
