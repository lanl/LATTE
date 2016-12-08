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

void M_InitWithLocal(Matrix &A, REAL *Local, int M, int N) {

  int Mpad, Npad;
  int cdev;
  int PadSize;
  // Get current device
  cudaGetDevice(&cdev);

  A.M=M;
  A.N=N;

  PadSize = nblocks*16;

//  if ( nblocks <= 2 ) {	
    
    if ( sizeof(REAL) == 8) {
      
      // double precision case - multiples of 64 best
      
      if (N <= 736) {
	
	Npad = nblocks*((N-1)/nblocks + 1);
	
      } else if ( N > 736 ) {
	
	Npad = PadSize*((N-1)/PadSize + 1);
	
      } 
      
    }
    
    if (sizeof(REAL) == 4) {
      
      // Single precision dimensions 
      
      if (N <= 448) {
	
        Npad = nblocks*((N-1)/nblocks + 1);
	
      } else if ( N > 448 ) {
	
        Npad = PadSize*((N-1)/PadSize + 1);
	
      }
      
      
    }   
      
       

//  } else if ( nblocks > 2) {

    
//    Npad = nblocks*((N - 1)/nblocks + 1);

//    printf("%d %d %d\n", nblocks, N, Npad);

//  } 
  
//   Npad = N;
  
  A.DN = Npad;
  A.DM = A.DN*M/N;

//  printf("InitWithLocal: %d %d %d %d \n", M, Mpad, A.DM, A.DN);  

  
  A.Local=Local;

  for (int d = 0; d < ndevices; d++) {
    cudaSetDevice(d);
    cudaMalloc((void **)&A.Device[d], A.DM*A.DN*sizeof(REAL));
    cudaMemset(A.Device[d], '\0', A.DM*A.DN*sizeof(REAL));
  }

  // Restore device
  cudaSetDevice(cdev);
}
