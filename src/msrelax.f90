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

SUBROUTINE MSRELAX

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE COULOMBARRAY
  USE VIRIALARRAY
  USE RELAXCOMMON
  USE MYPRECISION
  USE CONSTRAINTS_MOD

  IMPLICIT NONE

  INTEGER :: ITER, I, BREAKLOOP
  REAL(LATTEPREC) :: MAXF, DELTAF, PREVE, DELTAENERGY, PHI
  REAL(LATTEPREC) :: PSCALE
  REAL(LATTEPREC), PARAMETER  :: MAXPSCALE = 0.01
  IF (EXISTERROR) RETURN

  !  REAL(LATTEPREC) :: PIE = 3.141592654

  ITER = 0
  MAXF = 100000000000.0
  ENTE = ZERO


  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="monitorrelax.xyz")

  WRITE(20, 10) NATS
  WRITE(20, '("Molecular statics relaxation")')
  DO I = 1, NATS
     WRITE(20,11) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
  ENDDO

10 FORMAT(I6)
11 FORMAT(A2, 1X, 3F24.9)

  IF (CONTROL .EQ. 1) THEN
     CALL ALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL ALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIALLOCATE
  ENDIF

  IF (ELECTRO .EQ. 1) THEN
     CALL ALLOCATECOULOMB
     CALL INITCOULOMB
  ENDIF

  !
  ! Allocate stuff for building the neighbor lists
  !

  CALL ALLOCATENEBARRAYS

  CALL NEBLISTS(0)

  IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

  !  CALL BLDNEWHS(0)

  IF (KON .EQ. 0) THEN

     IF (SPONLY .EQ. 0) THEN
        CALL BLDNEWHS_SP
     ELSE
        CALL BLDNEWHS
     ENDIF

  ELSE

     CALL KBLDNEWH

  ENDIF


  IF (SPINON .EQ. 1) THEN
     CALL GETDELTASPIN
     CALL BLDSPINH
  ENDIF

  IF (CONTROL .EQ. 5) THEN
     CALL GERSHGORIN
     CALL SP2FERMIINIT
  ENDIF

  IF ( RELTYPE .EQ. "CG" ) ALLOCATE(OLDF(3,NATS), D1(3,NATS), OLDD(3,NATS))


  WRITE(6,'("# Iteration       Max. Force         Total Energy        Pressure")')

  BREAKLOOP  = 0

  IF (RELTYPE .EQ. "SD" .OR. RELTYPE .EQ. "CG") THEN

     DO WHILE ( BREAKLOOP .EQ. 0)

        ITER = ITER + 1

        IF (ELECTRO .EQ. 0) THEN
           CALL QNEUTRAL(0, 1)  ! Local charge neutrality
        ELSE
           CALL QCONSISTENCY(0,1) ! Self-consistent charges
        ENDIF

        IF (ELECTRO .EQ. 1) CALL GETCOULE

        ESPIN = ZERO
        IF (SPINON .EQ. 1) CALL GETSPINE

        IF (CONTROL .NE. 1 .AND. CONTROL .NE. 2 &
             .AND. KBT .GT. 0.000001 ) CALL ENTROPY

        CALL GETFORCE ! Get all the forces

        IF (FREEZE .EQ. 1) CALL FREEZE_ATOMS(FTOT)

        PREVE = TOTE

        CALL TOTENG

        IF (ELECTRO .EQ. 0) THEN

           TOTE = TRRHOH + EREP - ENTE

        ELSEIF (ELECTRO .EQ. 1) THEN

           TOTE = TRRHOH + EREP - ECOUL - ENTE

        ENDIF

        IF (SPINON .EQ. 1) TOTE = TOTE + ESPIN

        PREVF = MAXF
        CALL GETMAXF(MAXF)

        CALL GETPRESSURE


        WRITE(6,20) ITER, MAXF, TOTE, PRESSURE

        FLUSH(6)

20      FORMAT(1X,I6,10X,G14.6,1X,F20.10,1X,F6.2,1X,F12.6)

        DELTAENERGY = TOTE - PREVE
        DELTAF = MAXF - PREVF

        IF ( RELTYPE .EQ. "SD" ) THEN
           CALL STDESCENT(ITER, DELTAENERGY, DELTAF)
        ELSEIF ( RELTYPE .EQ. "CG" ) THEN
           CALL CONJGRADIENT(ITER, DELTAENERGY)
        ENDIF

        IF (ABS(MAXF) .LE. RLXFTOL .OR. ITER .GT. MXRLX ) BREAKLOOP = 1


        !
        ! After moving atoms, apply PBCs again
        !

        IF (MOD(ITER,25) .EQ. 0) CALL NEBLISTS(1)


        WRITE(20, 10) NATS
        WRITE(20, '("Molecular statics relaxation")')
        DO I = 1, NATS
           WRITE(20,11) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
        ENDDO

        IF (KON .EQ. 0) THEN

           IF (SPONLY .EQ. 0) THEN
              CALL BLDNEWHS_SP
           ELSE
              CALL BLDNEWHS
           ENDIF

        ELSE

           CALL KBLDNEWH

        ENDIF

     ENDDO

  ELSEIF (RELTYPE .EQ. "VL") THEN

     DO WHILE ( BREAKLOOP .EQ. 0)

        ITER = 0

        DO WHILE (ABS(MAXF) .GT. RLXFTOL .AND. ITER .LT. MXRLX)

           ITER = ITER + 1

           IF (ELECTRO .EQ. 0) THEN
              CALL QNEUTRAL(0,1)
           ELSE
              CALL QCONSISTENCY(0, 1)
           ENDIF

           IF (ELECTRO .EQ. 1) CALL GETCOULE

           ESPIN = ZERO
           IF (SPINON .EQ. 1) CALL GETSPINE

           IF ( CONTROL .NE. 1 .AND. CONTROL .NE. 2 &
                .AND. KBT .GT. 0.000001 ) CALL ENTROPY

           CALL GETFORCE

           IF (FREEZE .EQ. 1) CALL FREEZE_ATOMS(FTOT)

           CALL TOTENG

           PREVE = TOTE

           IF (ELECTRO .EQ. 0) THEN

              TOTE = TRRHOH + EREP - ENTE


           ELSEIF (ELECTRO .EQ. 1) THEN

              TOTE = TRRHOH + EREP - ECOUL - ENTE


           ENDIF

           IF (SPINON .EQ. 1) TOTE = TOTE + ESPIN

           PREVF = MAXF
           CALL GETMAXF(MAXF)

           CALL GETPRESSURE


           WRITE(6,21) ITER, MAXF, TOTE, PRESSURE, BOXDIMS(1)*BOXDIMS(2)*BOXDIMS(3)

21         FORMAT(1X,I6,10X,G14.6,1X,F20.10,1X,F6.2,1X,F12.6,1X,F12.6 )

           DELTAENERGY = TOTE - PREVE
           DELTAF = MAXF - PREVF

           CALL STDESCENT(ITER, DELTAENERGY, DELTAF)

           IF (MOD(ITER,25) .EQ. 0) THEN
              CALL NEBLISTS(1)
           ENDIF

           WRITE(20, 10) NATS
           WRITE(20, '("Molecular statics relaxation")')
           DO I = 1, NATS
              WRITE(20,11) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
           ENDDO

           IF (KON .EQ. 0) THEN

              IF (SPONLY .EQ. 0) THEN
                 CALL BLDNEWHS_SP
              ELSE
                 CALL BLDNEWHS
              ENDIF

           ELSE

              CALL KBLDNEWH

           ENDIF

        ENDDO

        CALL GETPRESSURE

        PSCALE = MIN(ABS(PRESSURE), MAXPSCALE)

        PSCALE = SIGN(PSCALE, PRESSURE)

        CR = CR * (ONE + PSCALE)

        BOX = BOX * (ONE + PSCALE)

        BOXDIMS = BOXDIMS * (ONE + PSCALE)

        CALL NEBLISTS(1)

        IF (KON .EQ. 0) THEN

           IF (SPONLY .EQ. 0) THEN
              CALL BLDNEWHS_SP
           ELSE
              CALL BLDNEWHS
           ENDIF

        ELSE

           CALL KBLDNEWH

        ENDIF

        IF (ABS(PRESSURE) .LT. 0.01) BREAKLOOP = 1

     ENDDO

  ENDIF

  CLOSE(20)

  CALL WRTRESTART(ITER)

  CALL WRTCFGS(-999)

  IF (CONTROL .EQ. 1) THEN
     !     CALL DEALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL DEALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIDEALLOCATE
  ENDIF

  CALL SUMMARY
  CALL FITTINGOUTPUT(0)

  IF (ELECTRO .EQ. 1)  CALL DEALLOCATECOULOMB

  IF (BASISTYPE .EQ. "NONORTHO") CALL DEALLOCATENONO

  CALL DEALLOCATENEBARRAYS

  IF ( RELTYPE .EQ. "CG" ) DEALLOCATE( OLDF, OLDD, D1 )


  RETURN

END SUBROUTINE MSRELAX
