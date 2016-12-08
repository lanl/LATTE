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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "Matrix.h"

extern cublasHandle_t handle;

void TestMultiply() {

  Matrix a, b;
  REAL trx, trx2;
  int hdim;
/*
  REAL *XZERO, *XONE;
  REAL YZERO = (REAL)0.0;
  REAL YONE = (REAL)1.0;
*/

  // Dim size
  hdim = 1000;

  // Init matrix
  M_Init(a, hdim, hdim);
  M_Init(b, hdim, hdim);

  // Send to GPU
  M_Push(a);
  M_Push(b);

  // Add identity matrix
  M_AddIdentity(a);
  M_AddIdentity(b);

  // Calculate trace
  trx = M_Trace(a);

  fprintf(stdout, " Test Trace a = %f \n", trx);

  // Calculate trace
  trx = M_Trace(b);

  fprintf(stdout, " Test Trace b = %f \n", trx);

  // Multiply by 2
  M_MultiplyScalar(&TWO,a);

#if REALSIZE==4
  //cublasSscal(handle, a.DM*a.DN, &TWO, a.Device, 1);
#elif REALSIZE==8
  //cublasDscal(handle, a.DM*a.DN, &TWO, a.Device, 1);
#endif

  // Calculate trace
  trx = M_Trace(a);

  fprintf(stdout, " Test MultiplyScalar Trace2 = %f \n", trx);

  // Multiply b^2
  //M_Multiply(b, b, a);
  M_Multiply(&ONE, b, b, &ZERO, a);
/*
  XZERO = &YZERO;
  XONE = &YONE;
#if REALSIZE==4
  cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, b.DM, b.DN, b.DN, XONE, b.Device, b.DM, b.Device, b.DM, XZERO, a.Device, a.DM);
#elif REALSI&ZE==8
  cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, b.DM, b.DN, b.DN, XONE, b.Device, b.DM, b.Device, b.DM, XZERO, a.Device, a.DM);
#endif
*/

  // Calculate trace
  trx = M_Trace(a);

  fprintf(stdout, " Test Multiply Trace4 = %f \n", trx);

  // Multiply
  M_Multiply(&ONE, b, b, &ZERO, a);

  // Calculate trace
  trx = M_Trace(a);

  fprintf(stdout, " Test Multiply Trace3 = %f \n", trx);

  // Multiply b^2
  M_Multiply(b, b, a);

  // Calculate trace
  trx = M_Trace(a);

  fprintf(stdout, " Test Multiply Trace4 = %f \n", trx);

  // MultiplyScalarSum
  M_MultiplyScalarSum(&TWO, b, a);

  trx = M_Trace(a);

  fprintf(stdout, " Test MultiplyScalarSum Trace5 = %f \n", trx);
 
  // Finish
  M_DeallocateLocal(a);
  M_DeallocateDevice(a);
  M_DeallocateLocal(b);
  M_DeallocateDevice(b);

}

/*
int main() {

  TestMultiply();

}
*/
