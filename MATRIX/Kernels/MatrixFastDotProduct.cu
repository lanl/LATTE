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

#include "Kernels.h"

__global__ void MatrixFastTraceX2Kernel(REAL *A, REAL *trace, int M, int num_threads) {
  
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  int i=idx, j=idx, active_threads, sum_index, num_elements=M*M, m, n;
  __shared__ REAL intermediate_sums[2][512];
  intermediate_sums[0][i]=0.0;
  while(j<num_elements) {
    m=j/M;
    n=j%M;
    intermediate_sums[0][i]+=(A[m+n*M]*A[n+m*M]);
    j+=num_threads;
  }
  active_threads=num_threads;
  sum_index=0;

  __syncthreads();

  while(active_threads>2) {
    if (i==0) intermediate_sums[sum_index][active_threads]=0.0;
    if (i==1) intermediate_sums[sum_index][active_threads+1]=0.0;
    __syncthreads();
    active_threads=active_threads/2+1;

    j=idx;
    if (j<active_threads) {
      intermediate_sums[!sum_index][j]=intermediate_sums[sum_index][j*2]+intermediate_sums[sum_index][j*2+1];
    }
    __syncthreads();
    sum_index=!sum_index;
  }
  if (i==0) {
    *trace=intermediate_sums[sum_index][0]+intermediate_sums[sum_index][1];
  }
}

