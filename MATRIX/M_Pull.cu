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

extern cublasHandle_t* handle;
extern int ndevices;
extern int nblocks;
extern cudaStream_t stream[];
extern cudaEvent_t event[];

void M_Pull(Matrix A) {

  cudaSetDevice(0);
  cublasGetMatrix(A.M, A.N, sizeof(REAL), A.Device[0], A.DM, A.Local, A.M);

}

void M_PullMgpu(Matrix A, int idevice) {

  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

  cudaSetDevice(idevice);
  cublasGetMatrix(A.M, A.N, sizeof(REAL), A.Device[idevice], A.DM, A.Local, A.M);

  // Restore device
  cudaSetDevice(cdev);
}

void M_PullMgpu(Matrix A) {

  int sub = A.DN / nblocks;
  int size = A.DN * sub;
  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

  cudaSetDevice(0);

  // For all GPUs
  if (ndevices > 1) {
    for (int d = 1; d < ndevices; ++d) {

      for (int b = d; b < nblocks; b += ndevices) {

        cudaMemcpy(A.Device[0]+b*size, A.Device[d]+b*size, size*sizeof(REAL), cudaMemcpyDefault);
      }
    }
  }

  // Move to CPU
  cublasGetMatrix(A.M, A.N, sizeof(REAL), A.Device[0], A.DM, A.Local, A.M);

  // Restore device
  cudaSetDevice(cdev);
}

void M_CollectDistributeMgpu(Matrix A) {

  int sub = A.DN / nblocks;
  int size = A.DN * sub;
  int cdev;

//  printf("%d %d %d %d \n", A.DN, A.DM, sub, size);

  // Get current device
  cudaGetDevice(&cdev);

//  cudaSetDevice(0);


  if (ndevices > 1) {

  for (int d = 0; d < ndevices; ++d) {

    for (int dp = 0; dp < ndevices; ++dp) {

       if (d != dp) {

          // Copy d's data to device dp 
	  
	  for (int b = d; b < nblocks; b += ndevices) {

	  cudaSetDevice(dp);
	  cudaMemcpy(A.Device[dp] + b*size, A.Device[d] + b*size, size*sizeof(REAL), cudaMemcpyDefault);

	  }

	}
      }

   }

   }

  // Collect blocks to GPU0 from other GPUs
/*  if (ndevices > 1) {
    for (int d = 1; d < ndevices; ++d) {

      for (int b = d; b < nblocks; b += ndevices) {

        cudaMemcpy(A.Device[0]+b*size, A.Device[d]+b*size, size*sizeof(REAL), cudaMemcpyDefault);
      }
    }
  }

  // Distribute Matrix to other GPUs
  if (ndevices > 1) {
    for (int d = 1; d < ndevices; ++d) {
      cudaSetDevice(d);
      cudaMemcpy(A.Device[d], A.Device[0], A.DM*A.DN*sizeof(REAL), cudaMemcpyDefault);
    }
  } */

  // Restore device
  cudaSetDevice(cdev);
}
