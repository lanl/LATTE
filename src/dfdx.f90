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

SUBROUTINE DFDX(I, J, L1, L2, M1, M2, R, ALPHA, COSBETA, WHICHINT, DFDA, DFDB, DFDR)

  ! Build derivative defined in PRB 72 165107 eq. (13)
  
  USE WIGNERARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, L1, L2, M1, M2, MP
  REAL(LATTEPREC) :: ALPHA, COSBETA, R, DFDA, DFDB, DFDR
  REAL(LATTEPREC), EXTERNAL :: WIGNERD, AM, BM
  REAL(LATTEPREC) :: WIG1, WIG2, SL1, SL2, AM1, AM2, BM1, BM2, INTEGRAL, TL1, TL2
  REAL(LATTEPREC) :: DSLA1, DSLB1, DSLA2, DSLB2, DTLA1, DTLB1, DTLA2, DTLB2 
  REAL(LATTEPREC) :: DWIG1, DWIGNEG1, DWIG2, DWIGNEG2
  REAL(LATTEPREC) :: WIGLIST1(-4:4), WIGLIST2(-4:4)
  REAL(LATTEPREC) :: BONDINT, DBONDINT
  CHARACTER(LEN=1) :: WHICHINT

  DO MP = -MIN(L1, L2)-1, MIN(L1, L2)+1
     WIGLIST1(MP) = WIGNERD(L1, ABS(M1), MP, COSBETA)
     WIGLIST2(MP) = WIGNERD(L2, ABS(M2), MP, COSBETA)
  ENDDO

  WIG1 = WIGLIST1(0)
  WIG2 = WIGLIST2(0)
  AM1 = AM(M1, ALPHA)
  AM2 = AM(M2, ALPHA)
  BM1 = BM(M1, ALPHA)
  BM2 = BM(M2, ALPHA)

  DWIG1 = HALF * MYSQRT(L1) * MYSQRT(L1 + 1) * &
       WIGLIST1(-1) - HALF * MYSQRT(L1) * MYSQRT(L1 + 1) * &
       WIGLIST1(1)

  DWIG2 = HALF * MYSQRT(L2) * MYSQRT(L2 + 1) * &
       WIGLIST2(-1) - HALF * MYSQRT(L2) * MYSQRT(L2 + 1) * &
       WIGLIST2(1)

  CALL DUNIVSCALE(I, J, L1, L2, 0, R, WHICHINT, BONDINT, DBONDINT)

  DFDA = TWO * WIG1 * WIG2 * BONDINT * &
       (ABS(M1) * BM1 * AM2 + ABS(M2) * &
       AM1 * BM2)
  
  DFDB = TWO * AM1 * AM2 * BONDINT * &
       (DWIG1 * WIG2 + WIG1 * DWIG2)

  DFDR = TWO * AM1 * AM2 * WIG1 * WIG2 * DBONDINT

  DO MP = 1, MIN(L1, L2)

     CALL DUNIVSCALE(I, J, L1, L2, MP, R, WHICHINT, BONDINT, DBONDINT)

     SL1 = AM1*(REAL(MINUSONEPOW(MP)) * WIGLIST1(MP) + WIGLIST1(-MP))
     SL2 = AM2*(REAL(MINUSONEPOW(MP)) * WIGLIST2(MP) + WIGLIST2(-MP))

     IF (M1 .EQ. 0) THEN
        TL1 = ZERO
        DTLA1 = ZERO
     ELSE 
        TL1 = BM1*(REAL(MINUSONEPOW(MP)) * WIGLIST1(MP) - WIGLIST1(-MP))
        DTLA1 = -ABS(M1)*AM1*(REAL(MINUSONEPOW(MP)) * WIGLIST1(MP) - WIGLIST1(-MP))
     ENDIF

     IF (M2 .EQ. 0) THEN
        TL2 = ZERO
        DTLA2 = ZERO
     ELSE 
        TL2 = BM2*(REAL(MINUSONEPOW(MP)) * WIGLIST2(MP) - WIGLIST2(-MP))
        DTLA2 = -ABS(M2)*AM2*(REAL(MINUSONEPOW(MP)) * WIGLIST2(MP) - WIGLIST2(-MP))
     ENDIF

     DSLA1 = ABS(M1)*BM1*(REAL(MINUSONEPOW(MP)) * WIGLIST1(MP) + WIGLIST1(-MP))
     DSLA2 = ABS(M2)*BM2*(REAL(MINUSONEPOW(MP)) * WIGLIST2(MP) + WIGLIST2(-MP))

     DWIG1 = HALF * MYSQRT(L1+MP) * MYSQRT(L1 - MP + 1) * &
          WIGLIST1(MP-1) - HALF * MYSQRT(L1 - MP) * MYSQRT(L1 + MP + 1) * &
          WIGLIST1(MP+1)

     DWIGNEG1 = HALF * MYSQRT(L1+(-MP)) * MYSQRT(L1 - (-MP) + 1) * &
          WIGLIST1((-MP)-1) - &
          HALF * MYSQRT(L1 - (-MP)) * MYSQRT(L1 + (-MP) + 1) * &
          WIGLIST1((-MP)+1)

     DWIG2 = HALF * MYSQRT(L2+MP) * MYSQRT(L2 - MP + 1) * &
          WIGLIST2(MP-1) - HALF * MYSQRT(L2 - MP) * MYSQRT(L2 + MP + 1) * &
          WIGLIST2(MP+1)

     DWIGNEG2 = HALF * MYSQRT(L2+(-MP)) * MYSQRT(L2 - (-MP) + 1) * &
          WIGLIST2((-MP)-1) - &
          HALF * MYSQRT(L2 - (-MP)) * MYSQRT(L2 + (-MP) + 1) * &
          WIGLIST2((-MP)+1)

     DSLB1 = AM1*(REAL(MINUSONEPOW(MP)) * DWIG1 + DWIGNEG1)
     DSLB2 = AM2*(REAL(MINUSONEPOW(MP)) * DWIG2 + DWIGNEG2)


     IF (M1 .EQ. 0) THEN
        DTLB1 = ZERO
     ELSE
        DTLB1 = BM1*(REAL(MINUSONEPOW(MP)) * DWIG1 - DWIGNEG1)
     ENDIF

     IF (M2 .EQ. 0) THEN
        DTLB2 = ZERO
     ELSE
        DTLB2 = BM2*(REAL(MINUSONEPOW(MP)) * DWIG2 - DWIGNEG2)
     ENDIF

     DFDA = DFDA + (DSLA1*SL2 + SL1*DSLA2 + DTLA1*TL2 + TL1*DTLA2)*BONDINT

     DFDB = DFDB + (DSLB1*SL2 + SL1*DSLB2 + DTLB1*TL2 + TL1*DTLB2)*BONDINT

     DFDR = DFDR + (SL1 * SL2 + TL1 * TL2) * DBONDINT

  ENDDO

  RETURN 

END SUBROUTINE DFDX
