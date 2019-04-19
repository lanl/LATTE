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

SUBROUTINE BOEVECSPRG

#ifdef PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE MYPRECISION
  USE MDARRAY
  USE NONOARRAY
  USE SPARSEARRAY

  USE BML
  USE PRG_DENSITYMATRIX_MOD

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: ITER, BREAKLOOP, LOOPTARGET
  REAL(LATTEPREC) :: OCCTARGET, OCC, FDIRAC, DFDIRAC
  REAL(LATTEPREC) :: OCCERROR, SHIFTCP, FDIRACARG, EXPARG
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: EBAND, QMIXORIG
  REAL(LATTEPREC) :: S, OCCLOGOCC_ELECTRONS, OCCLOGOCC_HOLES
  TYPE(BML_MATRIX_T) :: ORTHOH_BML, BO_BML

  IF (VERBOSE .GE. 1) WRITE(*,*) "In BOEVECSPRG ..."


  BO = ZERO

  ! OCCTARGET = BNDFIL*REAL(HDIM)

  !! Convert Hamiltonian to bml format
  !! H should be in orthogonal form, ORTHOH
  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, ORTHOH_BML)
  CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
       LATTEPREC, HDIM, HDIM, BO_BML)
  CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
       ORTHOH, ORTHOH_BML, ZERO, HDIM)

  IF (KBT .GT. 0.000001) THEN  ! This bit is for a finite electronic temperature

     CALL PRG_BUILD_DENSITY_T(ORTHOH_BML,BO_BML,NUMTHRESH,BNDFIL,KBT,CHEMPOT,EVALS)

  ELSE ! This bit is for zero electronic temperature

     CALL PRG_BUILD_DENSITY_T0(ORTHOH_BML,BO_BML,NUMTHRESH,BNDFIL,EVALS)

  ENDIF

  CALL BML_EXPORT_TO_DENSE(BO_BML, BO)

  CALL BML_DEALLOCATE(BO_BML)
  CALL BML_DEALLOCATE(ORTHOH_BML)


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

  ! BO = TWO * BO

  OCC =  BNDFIL*FLOAT(HDIM)


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

      ENDIF

      ENTE = -KBT*S

        ! Compute the gap only when we have to...

  IF (MOD(INT(TOTNE),2) .EQ. 0) THEN
     EGAP = EVALS(INT(OCC) + 1) - EVALS(INT(OCC))
  ELSE
     EGAP = ZERO
  ENDIF


  RETURN

#endif

END SUBROUTINE BOEVECSPRG
