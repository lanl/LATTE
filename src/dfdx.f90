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

  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, L1, L2, M1, M2, MP
  REAL(LATTEPREC) :: ALPHA, COSBETA, R, DFDA, DFDB, DFDR
  REAL(LATTEPREC), EXTERNAL :: SLMMP, TLMMP, AM, BM, WIGNERD, UNIVSCALE, DUNIVSCALE
  REAL(LATTEPREC), EXTERNAL :: DSLMMPDA, DTLMMPDA, DSLMMPDB, DTLMMPDB, DWIGNERDDB
  REAL(LATTEPREC) :: WIG1, WIG2, SL1, SL2, AM1, AM2, INTEGRAL, TL1, TL2
  CHARACTER(LEN=1) :: WHICHINT

  WIG1 = WIGNERD(L1, ABS(M1), 0, COSBETA)
  WIG2 = WIGNERD(L2, ABS(M2), 0, COSBETA)
  AM1 = AM(M1, ALPHA)
  AM2 = AM(M2, ALPHA)
  INTEGRAL = UNIVSCALE(I, J, L1, L2, 0, R, WHICHINT)

  DFDA = TWO * WIG1 * WIG2 * INTEGRAL * &
       (ABS(M1) * BM(M1, ALPHA) * AM2 + ABS(M2) * &
       AM1 * BM(M2, ALPHA))
  
  DFDB = TWO * AM1 * AM2 * INTEGRAL * &
       (DWIGNERDDB(L1, ABS(M1), 0, COSBETA) * WIG2 + &
       WIG1 * DWIGNERDDB(L2, ABS(M2), 0, COSBETA))

  DFDR = TWO * AM1 * AM2 * WIG1 * WIG2 * &
       DUNIVSCALE(I, J, L1, L2, 0, R, WHICHINT) 


  DO MP = 1, MIN(L1, L2)
     
     INTEGRAL = UNIVSCALE(I, J, L1, L2, MP, R, WHICHINT)
     SL1 = SLMMP(L1, M1, MP, ALPHA, COSBETA)
     SL2 = SLMMP(L2, M2, MP, ALPHA, COSBETA)
     TL1 = TLMMP(L1, M1, MP, ALPHA, COSBETA)
     TL2 = TLMMP(L2, M2, MP, ALPHA, COSBETA)

     DFDA = DFDA + (DSLMMPDA(L1, M1, MP, ALPHA, COSBETA) * &
          SL2 + SL1 * DSLMMPDA(L2, M2, MP, ALPHA, COSBETA) + &
          DTLMMPDA(L1, M1, MP, ALPHA, COSBETA) * &
          TL2 + TL1 * DTLMMPDA(L2, M2, MP, ALPHA, COSBETA)) * &
          INTEGRAL

     DFDB = DFDB + (DSLMMPDB(L1, M1, MP, ALPHA, COSBETA) * &
          SL2 + SL1 * DSLMMPDB(L2, M2, MP, ALPHA, COSBETA) + &
          DTLMMPDB(L1, M1, MP, ALPHA, COSBETA) * TL2 + &
          TL1 * DTLMMPDB(L2, M2, MP, ALPHA, COSBETA)) * &
          INTEGRAL

     DFDR = DFDR + (SL1 * SL2 + TL1 * TL2) * &
          DUNIVSCALE(I, J, L1, L2, MP, R, WHICHINT)

  ENDDO

  RETURN 

END SUBROUTINE DFDX
