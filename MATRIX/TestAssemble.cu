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
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "Matrix.h"

void TestAssemble() {

  Matrix a;
  REAL trx, trx2;
  REAL *b; 
  int hdim, nblock, sub;

  // Dim size
  hdim = 64;
  nblock = 2;
  sub = hdim / nblock;

  // Allocate b and zero
  b = (REAL*)malloc(hdim*hdim*sizeof(REAL));
  memset(b, '\0', hdim*hdim*sizeof(REAL));
  
  // Set values in b
  for (int i = 0; i < hdim*hdim; i++) {
    b[i] = i;
  }

  // Print out b
  printf("\n b = \n");
  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
        printf("%7.3f ", b[i*hdim+j]);
    }
    printf("\n");
  }
  printf("\n");

  // Init matrix
  M_InitWithLocal(a, b, hdim, hdim);

  // Send to GPU
  M_Push(a);

  // Reassemble
  M_AssembleMgpu(a, 0, nblock, sub, stream, ngpu);;

  // Bring back a
  M_Pull(a);

  // Print out reassembled a 
  printf("\n a = \n");
  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
        printf("%7.3f ", a.Local[i*hdim+j]);
    }
    printf("\n");
  }
  printf("\n");


  // Finish
  free(b);
  M_DeallocateLocal(a);
  M_DeallocateDevice(a);

}

