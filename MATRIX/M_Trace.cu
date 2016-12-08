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

extern int ndevices;
extern int nblocks;
extern cudaStream_t stream[];
extern cudaEvent_t event[];

REAL M_Trace(Matrix A) {

  // Size is N/2
  int size = A.DM >> 1;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);
  int smemSize = NUM_THREADS * sizeof(REAL);
  REAL *device_trace;
  REAL *local_trace = (REAL*)malloc(blockCount * sizeof(REAL));

  REAL trace=ZERO;

/*
  printf("N/2 = %d\n", size);
  printf("blockCount = %d\n", blockCount);
  printf("smemSize = %d\n ", smemSize);
  printf("NUM_THREADS = %d\n", NUM_THREADS);  
  printf("sizeof(REAL) = %ld\n", sizeof(REAL));
  printf("DM = %d\n", A.DM);
*/
 
  cudaMalloc(&device_trace, blockCount * sizeof(REAL));

  cudaSetDevice(0);

  MatrixFastTraceKernel<<<blockCount,NUM_THREADS,smemSize>>>(A.DM, A.DM, A.Device[0], device_trace, 0);

  // Copy to local variable
  cudaThreadSynchronize();
  cudaMemcpy(local_trace, device_trace, blockCount * sizeof(REAL), cudaMemcpyDeviceToHost);
  cudaFree(device_trace);

  for (int i = 0; i < blockCount; i++) {
    trace += local_trace[i];
  }

  free(local_trace);

  return trace;
}

// This does not work correctly
REAL M_TraceMgpu(Matrix A, int idevice) {

  // Size is N/2
  int sub = A.DN / ndevices;
  int bsize = A.DN * sub;
  int size = sub >> 1;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);
  int smemSize = NUM_THREADS * sizeof(REAL);
  REAL *device_trace;
  REAL *local_trace = (REAL*)malloc(blockCount * sizeof(REAL));

  REAL trace=ZERO;
  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

/*
  printf("N/2 = %d\n", size);
  printf("blockCount = %d\n", blockCount);
  printf("smemSize = %d\n ", smemSize);
  printf("NUM_THREADS = %d\n", NUM_THREADS);  
  printf("sizeof(REAL) = %ld\n", sizeof(REAL));
  printf("DM = %d\n", A.DM);
*/

  cudaMalloc(&device_trace, blockCount * sizeof(REAL));

  cudaSetDevice(idevice);

  MatrixFastTraceKernel<<<blockCount,NUM_THREADS,smemSize>>>(A.DM, sub, A.Device[idevice]+idevice*bsize, device_trace, idevice*sub);

  // Copy to local variable
  cudaThreadSynchronize();
  cudaMemcpy(local_trace, device_trace, blockCount * sizeof(REAL), cudaMemcpyDeviceToHost);
  cudaFree(device_trace);

  for (int i = 0; i < blockCount; i++) {
    trace += local_trace[i];
  }

  free(local_trace);

  // Restore device
  cudaSetDevice(cdev);

  return trace;
}

REAL M_TraceMgpu(Matrix A) {

  // Size is N * block size
  int sub = A.DN / nblocks;
  int bsize = A.DN * sub;
  int size = sub >> 1;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);
  int smemSize = NUM_THREADS * sizeof(REAL);
  REAL *device_trace;
  REAL *local_trace = (REAL*)malloc(blockCount * sizeof(REAL));

  REAL trace=ZERO;
  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

/*
  printf("N/2 = %d\n", size);
  printf("blockCount = %d\n", blockCount);
  printf("smemSize = %d\n ", smemSize);
  printf("NUM_THREADS = %d\n", NUM_THREADS);
  printf("sizeof(REAL) = %ld\n", sizeof(REAL));
  printf("DM = %d\n", A.DM);
*/

  // For all GPUs
  for (int d = 0; d < ndevices; ++d) {

    cudaSetDevice(d);

    cudaMalloc(&device_trace, blockCount * sizeof(REAL));

    for (int b = d; b < nblocks; b+=ndevices) {

      MatrixFastTraceKernel<<<blockCount,NUM_THREADS,smemSize>>>(A.DM, sub, A.Device[d]+b*bsize, device_trace, b*sub);

      // Copy to local variable
      cudaThreadSynchronize();
      cudaMemcpy(local_trace, device_trace, blockCount * sizeof(REAL), cudaMemcpyDeviceToHost);

      for (int i = 0; i < blockCount; i++) {
        trace += local_trace[i];
      }
    }
    cudaFree(device_trace);
  }

  free(local_trace);

  // Restore device
  cudaSetDevice(cdev);

  return trace;
}

