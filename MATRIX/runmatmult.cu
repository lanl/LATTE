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

void runmatmult(int  hdim, REAL *x0_pointer, REAL *h_pointer) {

  int nomult = 10;  
  Matrix xtmp, x0;	
  struct timeval ptime;
  timeval t1, t2;
  double elapsedTime;

  M_Init(xtmp, hdim, hdim);
  M_InitWithLocal(x0, x0_pointer, hdim, hdim);
  
  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      x0_pointer[i*hdim + j] = h_pointer[i*hdim + j];
    }
  }
    
  // Copy x0 to all GPUs  
  M_PushMgpu( x0 );

  M_CopyMgpu(x0, xtmp);

/*  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start); */
  
  gettimeofday(&t1, NULL);
  
  for ( int i = 0; i < nomult; i++) {
  M_MultiplyMgpu( &ONE, x0, x0, &ONE, xtmp);
  }

  gettimeofday(&t2, NULL);

  elapsedTime = (t2.tv_sec - t1.tv_sec)*1000.0;
  elapsedTime += (t2.tv_usec - t1.tv_usec)/1000.0;
  

/*  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  
  printf("%i %f \n", hdim, milliseconds/nomult); */
  printf("%i %f \n", hdim, elapsedTime/double(nomult));
  
  M_CopyMgpu(xtmp, x0);
  	
  M_PullMgpu(x0);
  
  M_DeallocateLocal(xtmp);
  M_DeallocateDevice(xtmp);
  M_DeallocateDevice(x0);

}
