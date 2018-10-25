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

SUBROUTINE GETRHO(MDITER)

  USE CONSTANTS_MOD

#ifdef PROGRESSON
  USE SP2PROGRESS
#endif


  IMPLICIT NONE

  INTEGER :: MDITER
  IF (EXISTERROR) RETURN

  ! This subroutine selects and calls the subroutines used to
  ! compute the density matrix based on the value of CONTROL

  ! CONTROL = 1 : DIAGONALIZATION
  ! CONTROL = 2 : SP2
  ! CONTROL = 3 : FERMI OPERATOR EXPANSION (FINITE T)
  ! CONTROL = 4, 5 : EXPERIMENTAL FINITE T SP2

  IF (CONTROL .EQ. 1) THEN
#ifdef PROGRESSON
    IF (SPINON .EQ. 0) THEN
      CALL BOEVECSPRG()
    ELSE
      CALL DIAGMYH()
      CALL SPINRHOEVECS
      IF(VERBOSE >= 0) WRITE(*,*)"This is using the original LATTE routine. Spin-polarized non yet with PROGRESS/BML"
    ENDIF
#else
     CALL DIAGMYH()

     IF (SPINON .EQ. 0) THEN
        CALL BOEVECS()
     ELSE
        CALL SPINRHOEVECS
     ENDIF
#endif

  ELSEIF (CONTROL .EQ. 2) THEN

     IF  (SPARSEON .EQ. 0) THEN

        CALL GERSHGORIN

        CALL SP2PURE

        !        IF (MDITER .LE. 10) THEN
        !           CALL SP2GAP_SETUP
        !        ELSE
        !           CALL SP2GAP
        !        ENDIF

     ELSEIF (SPARSEON .EQ. 1) THEN

        !! CALL SP2PURE_SPARSE

#ifdef PROGRESSON

        CALL SP2PRG

#else

        IF (MDITER .LE. 10) THEN
           CALL SP2PURE_SPARSE_PARALLEL(MDITER)
        ELSE
           CALL SP2PURE_SPARSE_PARALLEL_SIMPLE(MDITER)
        ENDIF

#endif

     ELSEIF (SPARSEON .EQ. 2) THEN

        IF (MDITER .LE. 10) THEN
           CALL SP2PURE_SPARSE_PARALLEL(MDITER)
        ELSE
           CALL SP2PURE_SUBGRAPH_PARALLEL(MDITER)
        ENDIF

     ENDIF


     !     ELSEIF (SPARSEON .EQ. 1) THEN

     !        CALL SP2PURE_SPARSE
     !        CALL SP2PURE_SPARSE_PARALLEL(MDITER)
     !     ELSEIF (SPARSEON .EQ. 2) THEN
     !        CALL SP2PURE_SPARSE_PARALLEL(MDITER)
     !     ENDIF

  ELSEIF (CONTROL .EQ. 3) THEN

     CALL FERMIEXPANS

  ELSEIF (CONTROL .EQ. 4) THEN

     CALL GERSHGORIN
     CALL SP2T

  ELSEIF (CONTROL .EQ. 5) THEN

     CALL SP2FERMI

  ENDIF

  RETURN

END SUBROUTINE GETRHO
