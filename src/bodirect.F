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

SUBROUTINE BOEVECS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE MYPRECISION
  USE MDARRAY

#ifdef PROGRESSON
  USE PRG_DENSITYMATRIX_MOD
#endif

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: ITER, BREAKLOOP, LOOPTARGET
  REAL(LATTEPREC) :: OCCTARGET, OCC, FDIRAC, DFDIRAC
  REAL(LATTEPREC) :: OCCERROR, SHIFTCP, FDIRACARG, EXPARG
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: EBAND, QMIXORIG
  REAL(LATTEPREC) :: S, OCCLOGOCC_ELECTRONS, OCCLOGOCC_HOLES

  IF (EXISTERROR) RETURN

  BO = ZERO

  OCCTARGET = BNDFIL*REAL(HDIM)

  !  PRINT*, TOTNE, OCCTARGET

  ITER = 0

  BREAKLOOP = 0

  OCCERROR = 1000000000.0

  !
  ! The do-while loop uses a Newton-Raphson optimization of the chemical
  ! potential to obtain the correct occupation
  !

  IF (KBT .GT. 0.000001) THEN  ! This bit is for a finite electronic temperature

     IF(VERBOSE >= 2)WRITE(*,*)"Total charge =",SUM(DELTAQ)

#ifdef PROGRESSON

     CALL PRG_GET_FLEVEL(EVALS,KBT,BNDFIL,BREAKTOL,CHEMPOT)

#else

     DO WHILE (ABS(OCCERROR) .GT. BREAKTOL .AND. ITER .LT. 100)

        ITER = ITER + 1
        OCC = ZERO
        DFDIRAC = ZERO

        DO I = 1, HDIM

           FDIRACARG = (EVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           EXPARG = EXP(FDIRACARG)
           FDIRAC = ONE/(ONE + EXPARG)
           OCC = OCC + FDIRAC
           DFDIRAC = DFDIRAC + EXPARG*FDIRAC*FDIRAC

        ENDDO

        DFDIRAC = DFDIRAC/KBT

        OCCERROR = OCCTARGET - OCC

        IF (ABS(DFDIRAC) .LT. NUMLIMIT) DFDIRAC = SIGN(NUMLIMIT, DFDIRAC)

        SHIFTCP = OCCERROR/DFDIRAC

        IF (ABS(SHIFTCP) .GT. MAXSHIFT) SHIFTCP = SIGN(MAXSHIFT, SHIFTCP)

        CHEMPOT = CHEMPOT + SHIFTCP

        !        PRINT*, CHEMPOT, OCCERROR

        IF(VERBOSE >= 2)WRITE(*,*)"Occupation error = ",OCCERROR," Chemical potential = ",CHEMPOT

     ENDDO

     IF (ITER .EQ. 100) THEN
        CALL ERRORS("bodirect","Newton-Raphson scheme to find the Chemical potential does not converge")
        IF (EXISTERROR) RETURN
     ENDIF

     ! Now we have the chemical potential we can build the density matrix
#endif

     S = ZERO

     IF (MDON .EQ. 0 .OR. &
          (MDON .EQ. 1 .AND. MOD(ENTROPYITER, WRTFREQ) .EQ. 0 )) THEN

        DO I = 1, HDIM

           FDIRACARG = (EVALS(I) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           FDIRAC = ONE/(ONE + EXP(FDIRACARG))

           OCCLOGOCC_ELECTRONS = FDIRAC * LOG(FDIRAC)
           OCCLOGOCC_HOLES = (ONE - FDIRAC) * LOG(ONE - FDIRAC)

           S = S + TWO*(OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES)

        ENDDO

        ! Compute the gap only when we have to...

        ! If we have an even number of electrons

        IF (MOD(INT(TOTNE),2) .EQ. 0) THEN
           EGAP = EVALS(INT(OCCTARGET) + 1) - EVALS(INT(OCCTARGET))
        ELSE
           EGAP = ZERO
        ENDIF

     ENDIF

     ENTE = -KBT*S

     DO I = 1, HDIM

        FDIRACARG = (EVALS(I) - CHEMPOT)/KBT

        FDIRACARG = MAX(FDIRACARG, -EXPTOL)
        FDIRACARG = MIN(FDIRACARG, EXPTOL)

        FDIRAC = ONE/(ONE + EXP(FDIRACARG))

#ifdef DOUBLEPREC
        CALL DGER(HDIM, HDIM, FDIRAC, EVECS(:,I), 1, EVECS(:,I), 1, BO, HDIM)
#elif defined(SINGLEPREC)
        CALL SGER(HDIM, HDIM, FDIRAC, EVECS(:,I), 1, EVECS(:,I), 1, BO, HDIM)
#endif

     ENDDO

  ELSE ! This bit is for zero electronic temperature

     IF (MOD(INT(TOTNE),2) .NE. 0) THEN
        CALL ERRORS("bodirect","Odd number of electrons - run a spin-polarized calculation &
             & or use a finite electron temperature")
     ENDIF

     !
     ! This definition of the chemical potential is a little arbitrary
     !

     LOOPTARGET = NINT(OCCTARGET)

     CHEMPOT = HALF*(EVALS(LOOPTARGET) + EVALS(LOOPTARGET + 1))
     EGAP = EVALS(LOOPTARGET + 1) - EVALS(LOOPTARGET)

     DO I = 1, LOOPTARGET

#ifdef DOUBLEPREC
        CALL DGER(HDIM, HDIM, ONE, EVECS(:,I), 1, EVECS(:,I), 1, BO, HDIM)
#elif defined(SINGLEPREC)
        CALL SGER(HDIM, HDIM, ONE, EVECS(:,I), 1, EVECS(:,I), 1, BO, HDIM)
#endif
     ENDDO

  ENDIF

  IF (MDON .EQ. 1 .AND. MDADAPT .EQ. 1) THEN

     FULLQCONV = 0

     IF (EGAP .LT. 1.0D0) THEN
        FULLQCONV = 1
        MDMIX = 0.1
     ELSE
        QITER = 1
        MDMIX = 0.25
     ENDIF

  ENDIF

  BO = TWO * BO

  RETURN

END SUBROUTINE BOEVECS
