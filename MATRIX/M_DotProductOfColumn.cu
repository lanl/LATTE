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

REAL M_DotProductOfColumn(int i, Matrix A, Matrix B) {

  REAL dot=ZERO;

  // M/2
  int size = A.M >> 1;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);
  int smemSize = NUM_THREADS * sizeof(REAL);
  REAL *device_dot;
  REAL *local_dot = (REAL*)malloc(blockCount * sizeof(REAL));
  
  cudaMalloc(&device_dot, blockCount * sizeof(REAL));

  cudaSetDevice(0);

  FastDotProductKernel<<<blockCount,NUM_THREADS,smemSize>>>(A.M, &A.Device[0][i*A.DM], &B.Device[0][i*B.DM], device_dot);

  // Copy to local variable
  cudaThreadSynchronize();
  cudaMemcpy(local_dot, device_dot, blockCount * sizeof(REAL), cudaMemcpyDeviceToHost);
  cudaFree(device_dot);

  for (int i = 0; i < blockCount; i++)
    dot += local_dot[i];

  free(local_dot);

  return dot;
}
