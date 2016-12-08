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


void M_TransBlk2(int size, Matrix A, Matrix A2, int d) {

  int blockCount = (int) ceil((float)size/(float)NUM_THREADS);

//  A2.DM = A.DN;
//  A2.DN = A.DM;

  TransKernel<<<blockCount, NUM_THREADS,0, stream[d]>>>(size, A.DM, A.DN, A.Device[d], A2.Device[d]);

//  printf("%d %d \n", A2.DM, A2.DN);

//     if ( d == 0 ) {
/*      M_Pull(A);

        printf("\n A = \n");
  for (int i = 0; i < A.DM; i++) {
    for (int j = 0; j < A.DN; j++) {
        printf("%18.9f ", A.Local[i*A.N+j]);
    }
    printf("\n");
  }
  printf("\n");  

  M_Pull(A2);

        printf("\n T = \n");
  for (int i = 0; i < A2.DM; i++) {
    for (int j = 0; j < A2.DN; j++) {
        printf("%18.9f ", A2.Local[i*A2.N+j]);
    }
    printf("\n");
  }
  printf("\n");  
*/




//  } 


  




  

}



