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

FUNCTION DFDA(I, J, L1, L2, M1, M2, R, ALPHA, COSBETA, WHICHINT)

  ! Build derivative defined in PRB 72 165107 eq. (13)

  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, L1, L2, M1, M2, MP
  REAL(LATTEPREC) :: ALPHA, COSBETA, R, DFDA
  REAL(LATTEPREC), EXTERNAL :: SLMMP, TLMMP, AM, BM, WIGNERD, UNIVSCALE
  REAL(LATTEPREC), EXTERNAL :: DSLMMPDA, DTLMMPDA
  CHARACTER(LEN=1) :: WHICHINT

  DFDA = TWO * WIGNERD(L1, ABS(M1), 0, COSBETA) * &
       WIGNERD(L2, ABS(M2), 0, COSBETA) * &
       UNIVSCALE(I, J, L1, L2, 0, R, WHICHINT) * &
       (ABS(M1) * BM(M1, ALPHA) * AM(M2, ALPHA) + ABS(M2) * &
       AM(M1, ALPHA) * BM(M2, ALPHA))

  DO MP = 1, MIN(L1, L2)

     DFDA = DFDA + (DSLMMPDA(L1, M1, MP, ALPHA, COSBETA) * &
          SLMMP(L2, M2, MP, ALPHA, COSBETA) + &
          SLMMP(L1, M1, MP, ALPHA, COSBETA) * &
          DSLMMPDA(L2, M2, MP, ALPHA, COSBETA) + &
          DTLMMPDA(L1, M1, MP, ALPHA, COSBETA) * &
          TLMMP(L2, M2, MP, ALPHA, COSBETA) + &
          TLMMP(L1, M1, MP, ALPHA, COSBETA) * &
          DTLMMPDA(L2, M2, MP, ALPHA, COSBETA)) * &
          UNIVSCALE(I, J, L1, L2, MP, R, WHICHINT)

  ENDDO

  RETURN 

END FUNCTION DFDA
