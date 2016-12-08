/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

#include "Matrix.h"

// The Fortran interface

// Modified by MJC so we also send over FERMIM

extern "C" void solve_matrix_cg_(void *bo_ptr, void *rhoup_ptr, void *rhodown_ptr, int *hdim, void *cgtol2, int *spinflag, int *prec, int *fermim) {
  if (*prec==4) {
#if REALSIZE==4
    if (*spinflag) {
      solve_matrix_cg((float *)rhoup_ptr, *hdim, *((float *)cgtol2), *fermim);
      solve_matrix_cg((float *)rhodown_ptr, *hdim, *((float *)cgtol2), *fermim);
    }
    else {
      solve_matrix_cg((float *)bo_ptr, *hdim, *((float *)cgtol2), *fermim);
    }
#endif
  }
  if (*prec==8) {
#if REALSIZE==8
    if (*spinflag) {
      solve_matrix_cg((double *)rhoup_ptr, *hdim, *((double *)cgtol2), *fermim);
      solve_matrix_cg((double *)rhodown_ptr, *hdim, *((double *)cgtol2), *fermim);
    }
    else {
      solve_matrix_cg((double *)bo_ptr, *hdim, *((double *)cgtol2), *fermim);
    }
#endif
  }
}
