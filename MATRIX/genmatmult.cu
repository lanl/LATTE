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

#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>

#include "Matrix.h"

extern int ndevices;
extern int nblocks;

void genmatmult(int hdim, int tposea, int tposeb, REAL alpha, REAL beta, REAL *amat_pointer, REAL *bmat_pointer, REAL *cmat_pointer) {
  //void runmatmult(int  hdim, REAL *x0_pointer, REAL *h_pointer) {

  Matrix amat, bmat, cmat;

  M_InitWithLocal(amat, amat_pointer, hdim, hdim);
  M_InitWithLocal(bmat, bmat_pointer, hdim, hdim);
  M_InitWithLocal(cmat, cmat_pointer, hdim, hdim);

  // Copy Matrices to all GPUs. We only copy C if beta > 0  

  M_Push( amat );
  M_Push( bmat );	

  if (fabs(beta) > 1.0e-6) M_Push( cmat );
  
  M_Multiply(tposea, tposeb, &alpha, amat, bmat, &beta, cmat);
    	
  M_Pull(cmat);

  M_DeallocateDevice(amat);
  M_DeallocateDevice(bmat);
  M_DeallocateDevice(cmat);

}
