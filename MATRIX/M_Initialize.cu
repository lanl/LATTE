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
#include <stdio.h>

#include "Matrix.h"

// CUBLAS handle - required
cublasHandle_t *handle;
int maxdevices;
int ndevices;
int nblocks;
cudaStream_t stream[MAX_GPUS];
cudaEvent_t event[MAX_GPUS];

void M_Initialize( int NGPU ) {

  cublasStatus_t status;

  // Get number of GPUs
  cudaGetDeviceCount(&maxdevices);

  // Number of devices being used
  ndevices = NGPU;
//    fprintf(stdout, "Number of GPUs = %d \n", ndevices);

  // Number of blocks to spit up matrix into

  // Having the number of blocks = number of devices is usually a safe bet

    nblocks = 1*ndevices ;


  // Allow GPU 0 to access other GPUs
  if (ndevices > 1) {
    int accessible = 0;
    for (int d = 0; d < ndevices; d++) {
      for (int d2 = 0; d2 < ndevices; d2++) {
 //     printf("%d %d\n", d, d2);	
      
        if (d == d2) continue;

//	printf("2 %d %d\n", d, d2);	
        cudaDeviceCanAccessPeer(&accessible, d, d2);
//	printf("2 %d %d %d\n", d, d2, accessible);	
//        if (accessible == 0) {
	  if (accessible) {
          cudaSetDevice(d);
          cudaDeviceEnablePeerAccess(d2, 0);
//          printf("GPU %d can access GPU %d\n", d, d2);
        }
      }
    }      
  }

    // Create and assign steams and events
  for (int d = 0; d < ndevices; ++d) {
    cudaSetDevice(d);
    cudaStreamCreate(&stream[d]);
    cudaEventCreate(&event[d]);
  }
  cudaSetDevice(0); 

  // Initialize per GPU with it's own handle
  // Associate steam with handle
  handle = (cublasHandle_t*)malloc(ndevices * sizeof(cublasHandle_t));
  for (int d = 0; d < ndevices; ++d) {
    cudaSetDevice(d);

    status=cublasCreate(&(handle[d]));
    if (status!=CUBLAS_STATUS_SUCCESS) {
      if (status==CUBLAS_STATUS_ALLOC_FAILED) {
        printf("Could not allocate resources for GPU %d!\n", d);
      }
      printf("CuBLAS init failedfor GPU %d!\n", d);
      exit(-1);
    }

    // Associate stream with cublas handle
    cublasSetStream(handle[d], stream[d]);
  }

  // Create and assign steams and events
/*  for (int d = 0; d < ndevices; ++d) {
    cudaSetDevice(d);
    cudaStreamCreate(&stream[d]);
    cudaEventCreate(&event[d]);
  } */
 
  cudaSetDevice(0);

}

