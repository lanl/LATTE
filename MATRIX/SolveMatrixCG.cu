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

#include <stdio.h>

#include "Matrix.h"

void solve_matrix_cg(REAL *bo_ptr, int hdim, REAL cgtol2, int fermim) {
  int iter, breakloop, i;
  REAL error2;

  REAL xalpha, xbeta;
  REAL r0vec, p0vec, r1vec;

  Matrix bo, x2, a, tmpmat, r0, p0;

  M_InitWithLocal(bo, bo_ptr, hdim, hdim);
  M_Init(x2, hdim, hdim);
  M_Init(a, hdim, hdim);
  M_Init(tmpmat, hdim, hdim);
  M_Init(r0, hdim, hdim);
  M_Init(p0, hdim, hdim);
 
  // Copy bo to GPU0 and GPU1 or use from GPU0
  M_Push(bo);

  /* Added by MJC to minimize CPU/GPU communication - we don't pull
     the matrix back from the GPU until we've done all FERMIM loops */

  for (i = 0; i < fermim; i++) { 

    // A = (X^2 - X) * 2.0 - I
    // R0 = A * X - X^2
    // P0 = -R0
    M_Multiply( bo, bo, x2);
    M_Subtract(x2, bo, tmpmat);
    M_Multiply(TWO, tmpmat, a);
    M_AddIdentity(a);
    
    M_Multiply( a, bo, tmpmat);
    M_Subtract(tmpmat, x2, r0);
    M_Multiply(MINUS1, r0, p0);
    
    iter = 0;
    breakloop = 0;

    // Dot product of r0 and r0
    r0vec = M_DotProduct(r0, r0);

    while (breakloop == 0 && iter < 100) {

      iter++;

      // A * P0 - intermediate term used in CG
      M_Multiply( a, p0, tmpmat);

      // Dot product of P0 and A*P0
      p0vec = M_DotProduct( p0, tmpmat);

      // Set alpha
      if (p0vec > ZERO) xalpha = r0vec / p0vec;
      else xalpha = ZERO;

      // Do on GPU1
      // Copy p0 to GPU1 or use from GPU0
      // Or use bo and p) from GPU0
      // Calculate bo = alpha * p0 + bo
      M_MultiplyScalarSum( &xalpha, p0, bo);

      // Calculate r0 = alpha * tmpmat + r0
      M_MultiplyScalarSum( &xalpha, tmpmat, r0);

      // Dot product of r0 and r0
      r1vec = M_DotProduct( r0, r0);

      // Current error
      error2 = r1vec;

      // Set beta
      if (r0vec > ZERO) xbeta = r1vec / r0vec;
      else xbeta = ZERO;

      // Calculate p0 = beta * p0 - r0
      // p0 = beta * p0
      M_MultiplyScalar( &xbeta, p0);

      // p0 = -1.0 * r0 + p0
      M_MultiplyScalarSum( &MINUS1, r0, p0);

      //printf("iter = %d error2 = %e cgtol2= %e \n", iter, error2, cgtol2);

      if (error2 < cgtol2) breakloop = 1;

      r0vec = r1vec;
     
    }
    
  }

  if (iter == 100)
    printf("Solve Matrix CG reached 100 iterations: something is wrong! \n");

  //printf("iter = %d error2 = %e cgtol2= %e \n", iter, error2, cgtol2);

  M_Pull(bo);

  M_DeallocateDevice(bo);
  M_DeallocateLocal(x2);
  M_DeallocateDevice(x2);
  M_DeallocateLocal(a);
  M_DeallocateDevice(a);
  M_DeallocateLocal(tmpmat);
  M_DeallocateDevice(tmpmat);
  M_DeallocateLocal(p0);
  M_DeallocateDevice(p0);
  M_DeallocateLocal(r0);
  M_DeallocateDevice(r0);

}

