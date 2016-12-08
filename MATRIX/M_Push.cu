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

void M_Push(Matrix A) {

  cublasStatus_t status;

  cudaSetDevice(0);
  status=cublasSetMatrix(A.M, A.N, sizeof(REAL), A.Local, A.M, A.Device[0], A.DM);
  if (status!=CUBLAS_STATUS_SUCCESS) {
    if (status==CUBLAS_STATUS_NOT_INITIALIZED) printf("Push: CuBLAS not initialized!\n");
    if (status==CUBLAS_STATUS_INVALID_VALUE) printf("Push: Invalid value!\n");
    if (status==CUBLAS_STATUS_MAPPING_ERROR) printf("Push: error accessing GPU memory!\n");
    exit(1);
  }

}

void M_PushMgpu(Matrix A) {

  cublasStatus_t status;
  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

  for (int d = 0; d < ndevices; ++d) {

    cudaSetDevice(d);
    status=cublasSetMatrix(A.M, A.N, sizeof(REAL), A.Local, A.M, A.Device[d], A.DM);
    if (status!=CUBLAS_STATUS_SUCCESS) {
      if (status==CUBLAS_STATUS_NOT_INITIALIZED) printf("PushMgpu: CuBLAS not initialized!\n");
      if (status==CUBLAS_STATUS_INVALID_VALUE) printf("PushMgpu: Invalid value!\n");
      if (status==CUBLAS_STATUS_MAPPING_ERROR) printf("PushMgpu: error accessing GPU memory!\n");
      exit(1);
    }

  }

  // Restore device
  cudaSetDevice(cdev);

}

void M_PushDeviceMgpu(Matrix A) {

  int cdev;

  // Get current device
  cudaGetDevice(&cdev);

  for (int d = 1; d < ndevices; ++d) {

    cudaSetDevice(d);
    cudaMemcpy(A.Device[d], A.Device[0], A.DM*A.DN*sizeof(REAL), cudaMemcpyDefault);
  }

  // Restore device
  cudaSetDevice(cdev);
}

