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
extern "C" void sp2fermi_gpu_(void *bndfil, int  *hdim, int *spinon, void *bo_pointer, void *rhoup, void *rhodown, void *maxeval,
             void *h_pointer, void *hup, void *hdown, void *maxminusmin, void *chempot, int *norecs, void *signlist_pointer, void *beta_pointer, void *breaktol, int *prec) {
  if (*prec==4) {
#if REALSIZE==4
    if (*spinon) {
/*       sp2fermi_spin(*((float *)bndfil), *hdim, (float *)rhoup, (float *)rhodown,
                *((float *)maxeval), (float *)hup, (float *)hdown, *((float *)maxminusmin),
                (float *)chempot, *norecs, (float *)kbt_pointer, (float *)beta_pointer, *((float *)breaktol)); */
    }
    else {
        sp2fermi_nospin(*((float *)bndfil),     *hdim,              (float *)bo_pointer,
                        *((float *)maxeval),    (float *)h_pointer, *((float *)maxminusmin),
                        (float *)chempot,       *norecs,            (int *)signlist_pointer,
                        (float *)beta_pointer, *((float *)breaktol));
    }
#endif
  }
  if (*prec==8) {
#if REALSIZE==8
    if (*spinon) {
 /*       sp2fermi_spin(*((double *)bndfil), *hdim, (double *)rhoup, (double *)rhodown,
                *((double *)maxeval), (double *)hup, (double *)hdown, *((double *)maxminusmin),
                (double *)chempot, *norecs, (double *)kbt_pointer, (double *)beta_pointer, *((double *)breaktol)); */
    }
    else {
        sp2fermi_nospin(*((double *)bndfil),     *hdim,              (double *)bo_pointer,
                        *((double *)maxeval),    (double *)h_pointer, *((double *)maxminusmin),
                        (double *)chempot,       *norecs,            (int *)signlist_pointer,
                        (double *)beta_pointer, *((double *)breaktol)); 
    }
#endif
  }
}

/*
extern "C" void sp2fermi_init_(void *bndfil, int  *hdim, int *spinon, void *bo_pointer, void *rhoup, void *rhodown, void *maxeval,
             void *h_pointer, void *hup, void *hdown, void *maxminusmin, void *chempot, int *norecs, void *kbt_pointer, void *beta0_pointer, void *breaktol, int *prec) {
  if (*prec==4) {
    if (*spinon) {
        sp2fermi_init_spin(*((float *)bndfil), *hdim, (float *)rhoup, (float *)rhodown,
                *((float *)maxeval), (float *)hup, (float *)hdown, *((float *)maxminusmin),
                (float *)chempot, *norecs, (float *)kbt_pointer, (float *)beta0_pointer, *((float *)breaktol));
    }
    else {
        sp2fermi_init_nospin(*((float *)bndfil), *hdim, (float *)bo_pointer,
                *((float *)maxeval), (float *)h_pointer, *((float *)maxminusmin),
                (float *)chempot, *norecs, (float *)kbt_pointer, (float *)beta0_pointer, *((float *)breaktol));
    }
  }
  if (*prec==8) {
    if (*spinon) {
        sp2fermi_init_spin(*((double *)bndfil), *hdim, (double *)rhoup, (double *)rhodown,
                *((double *)maxeval), (double *)hup, (double *)hdown, *((double *)maxminusmin),
                (double *)chempot, *norecs, (double *)kbt_pointer, (double *)beta0_pointer, *((double *)breaktol));
    }
    else {
        sp2fermi_init_nospin(*((double *)bndfil), *hdim, (double *)bo_pointer,
                *((double *)maxeval), (double *)h_pointer, *((double *)maxminusmin),
                (double *)chempot, *norecs, (double *)kbt_pointer, (double *)beta0_pointer, *((double *)breaktol));
    }
  }
}
*/

