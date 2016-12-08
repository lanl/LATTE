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

// Assemble part of matrix for all GPUs
void M_AssembleMgpu(Matrix A, Matrix A2, int sub) {

  // Size is N * sub 
  int size = A.DN * sub;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);
  int cdev;

  // Save current device
  cudaGetDevice(&cdev);

/*
  printf("N = %d  sub = %d  nblocks = %d\n", A.DN, sub, nblocks);
  printf("size = %d  blockCount = %d  \n", size, blockCount);
*/
 
  // Do all blocks for each GPU
  for (int d = 0; d < ndevices; ++d) {

    cudaSetDevice(d);
    for (int b = d; b < nblocks; b += ndevices) {

      MatrixAssembleKernel<<<blockCount,NUM_THREADS,0,stream[d]>>>(size, A.DN, A.Device[d]+b*size, A2.Device[d]+b*size, nblocks, sub);

    }

  }

  // Wait till all reassembles are done
  M_Wait();

  // Reset device
  cudaSetDevice(cdev);

}

// Assemble part of a matrix for a single GPU
void M_AssembleMgpu(Matrix A, Matrix A2, int sub, int d) {

  // Size is N * sub
  int size = A.DN * sub;
  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);


//  printf("N = %d  sub = %d  nblocks = %d\n", A.DN, sub, nblocks);
//  printf("size = %d  blockCount = %d  \n", size, blockCount);

/*    if ( d == 0 ) {
      M_Pull(A2);

        printf("\n A2 = \n");
  for (int i = 0; i < A2.M; i++) {
    for (int j = 0; j < A2.N; j++) {
        printf("%18.9f ", A2.Local[i*A2.M+j]);
    }
    printf("\n");
  }
  printf("\n");  
  } */


  for (int b = d; b < nblocks; b += ndevices) {

//  printf("b = %d \n", b);

    MatrixAssembleKernel<<<blockCount,NUM_THREADS,0,stream[d]>>>(size, A.DN, A.Device[d]+b*size, A2.Device[d]+b*size, nblocks, sub);

  }

}



