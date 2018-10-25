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

SUBROUTINE INITIATEDM

  USE CONSTANTS_MOD
  USE SETUPARRAY
!  USE NEBLISTARRAY
  USE DMARRAY
!  USE NONOARRAY
!  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, CNT
!  INTEGER :: IBRA, IKET, LBRA, LKET, MBRA, MKET
!  INTEGER :: CNT, INDI, INDJ
!  INTEGER :: SWITCH, PREVJ
!  INTEGER :: PBCI, PBCJ, PBCK
!  INTEGER :: BASISI(5), BASISJ(5), LBRAINC, LKETINC
!  REAL(LATTEPREC) :: ALPHA, BETA, COSBETA, PHI, TMP, PERM
!  REAL(LATTEPREC) :: RIJ(3), MAGR2, MAGR, MAGRP, RCUTTB
!  REAL(LATTEPREC) :: MAXRCUT, MAXRCUT2
!  REAL(LATTEPREC) :: ANGFACTOR, AMMBRA, WIGLBRAMBRA
  !  REAL(LATTEPREC), ALLOCATABLE :: SAVEHDIAG(:)
!  REAL(LATTEPREC), EXTERNAL :: UNIVSCALE, WIGNERD, SLMMP, TLMMP, AM, BM
!  IF (EXISTERROR) RETURN


  H = ZERO

  CNT = 1

  ! Build diagonal elements
  DO I = 1, NATS

     K = ELEMPOINTER(I)

     SELECT CASE(BASIS(K))

     CASE("s")

        H_INDEX_START(I) = CNT
        CNT = CNT + 1
        H_INDEX_END(I) = CNT - 1

     CASE("p")

        H_INDEX_START(I) = CNT
        CNT = CNT + 3
        H_INDEX_END(I) = CNT - 1

     CASE("d")

        H_INDEX_START(I) = CNT
        CNT = CNT + 5
        H_INDEX_END(I) = CNT - 1

     CASE("f")

        H_INDEX_START(I) = CNT
        CNT = CNT + 7
        H_INDEX_END(I) = CNT - 1

     CASE("sp")

        H_INDEX_START(I) = CNT
        CNT = CNT + 4
        H_INDEX_END(I) = CNT - 1

     CASE("sd")

        H_INDEX_START(I) = CNT
        CNT = CNT + 6
        H_INDEX_END(I) = CNT - 1

     CASE("sf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 8
        H_INDEX_END(I) = CNT - 1

     CASE("pd")

        H_INDEX_START(I) = CNT
        CNT = CNT + 8
        H_INDEX_END(I) = CNT - 1

     CASE("pf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 10
        H_INDEX_END(I) = CNT - 1

     CASE("df")

        H_INDEX_START(I) = CNT
        CNT = CNT + 12
        H_INDEX_END(I) = CNT - 1

     CASE("spd")

        H_INDEX_START(I) = CNT
        CNT = CNT + 9
        H_INDEX_END(I) = CNT - 1

     CASE("spf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 11
        H_INDEX_END(I) = CNT - 1

     CASE("sdf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 13
        H_INDEX_END(I) = CNT - 1

     CASE("pdf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 15
        H_INDEX_END(I) = CNT - 1

     CASE("spdf")

        H_INDEX_START(I) = CNT
        CNT = CNT + 16
        H_INDEX_END(I) = CNT - 1

     END SELECT

  ENDDO

END SUBROUTINE INITIATEDM
