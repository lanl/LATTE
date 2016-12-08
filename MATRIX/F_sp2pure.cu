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

// Here is our fortran interface
extern "C" void sp2purify_(void *bndfil, int  *hdim, int *spinon, void *bo_pointer, void *rhoup, void *rhodown, void *maxeval,
             void *h_pointer, void *hup, void *hdown, void *maxminusmin, int *minsp2iter, int* sp2convint, int *prec) {
  if (*spinon) {
    if (*prec==4) {
#if REALSIZE==4
      sp2pure_spin3(*((float *)bndfil), *hdim, (float *)rhoup, (float *)rhodown,
       *((float *)maxeval), (float *)hup, (float *)hdown,
       *((float *)maxminusmin), *minsp2iter, *sp2convint);
#endif
    }
    if (*prec==8) {
#if REALSIZE==8
      sp2pure_spin3(*((double *)bndfil), *hdim, (double *)rhoup, (double *)rhodown,
       *((double *)maxeval), (double *)hup, (double *)hdown,
       *((double *)maxminusmin), *minsp2iter, *sp2convint);
#endif
    }
  }
  else {
    if (*prec==4) {
#if REALSIZE==4
      sp2pure_nospin3(*((float *)bndfil), *hdim, (float *)bo_pointer,
       *((float *)maxeval), (float *)h_pointer, *((float *)maxminusmin),
       *minsp2iter, *sp2convint);
#endif
    }
    if (*prec==8) {
#if REALSIZE==8
      sp2pure_nospin3(*((double *)bndfil), *hdim, (double *)bo_pointer,
       *((double *)maxeval), (double *)h_pointer, *((double *)maxminusmin),
       *minsp2iter, *sp2convint);
#endif
    }
  }
}


