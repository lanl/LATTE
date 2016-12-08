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

#include "Matrix.h"

extern int ndevices;
extern int nblocks;

void sp2pure_nospin3(REAL bndfil, int  hdim, REAL *x0_pointer, REAL maxeval,
             REAL *h_pointer, REAL maxminusmin, int minsp2iter, 
	     int sp2convint) {

  int iter, breakloop;
  REAL trx, trxtmp, trxold, trx0;
  REAL occ, limdiff;
  REAL idemperr = ZERO, idemperr1 = ZERO, idemperr2 = ZERO;
  REAL idemtol = SMALL_NUMBER;
  Matrix xtmp, x0;
  
  M_Init(xtmp, hdim, hdim);
  M_InitWithLocal(x0, x0_pointer, hdim, hdim);

  occ=bndfil*hdim;
  
  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      x0_pointer[i*hdim + j] = -h_pointer[i*hdim + j]/maxminusmin;
    }
  }
  
  for (int i = 0; i < hdim; i++) {
    x0_pointer[i*hdim + i] = maxeval/maxminusmin + x0_pointer[i*hdim + i];
  } 

  
  // Copy x0 to all GPUs  
  M_PushMgpu( x0 );

  trx =  M_TraceMgpu( x0 );

  //fprintf(stdout, "initial trace x0 = %f \n", trx);
  
  iter = 0;
  
  breakloop = 0;
  
  while( breakloop == 0 && iter < 100 ) {
    
    iter++;

    // On multiple GPUs: xtmp = x0 * x0 + x0

    M_CopyMgpu(x0, xtmp);

    M_MultiplyMgpu( &MINUS1, x0, x0, &ONE, xtmp);

      trxtmp = M_TraceMgpu( xtmp );

      trxold = trx;

//      fprintf(stdout, "%d trace xtmp = %f \n", iter, trxtmp);

#if REALSIZE == 4    
    limdiff = fabs(trx - trxtmp - occ) - fabs(trx + trxtmp - occ);
#elif REALSIZE == 8
    limdiff = abs(trx - trxtmp - occ) - abs(trx + trxtmp - occ);
#endif
    
    if (limdiff > idemtol ) {
      
      // x0 = x0 + xtmp
      M_MultiplyScalarSumMgpu( &ONE, xtmp, x0);
      trx = trx + trxtmp;
      
    } else if ( limdiff < -idemtol ) {
      
      // x0 = x0 - xtmp
      M_MultiplyScalarSumMgpu( &MINUS1, xtmp, x0);

      trx = trx - trxtmp;
      
    }

    // Distribute x0
    if (ndevices > 1) M_CollectDistributeMgpu(x0);

    idemperr2 = idemperr1;	
    idemperr1 = idemperr;
#if REALSIZE == 4
    idemperr = fabs( trx - trxold );
    //idemperr = fabs( trxtmp );
#elif REALSIZE == 8
    idemperr = abs( trx - trxold );
    //idemperr = abs( trxtmp );
#endif    
    
    if ( sp2convint == 0 && iter >= minsp2iter && 
	 (idemperr2 <= idemperr || idemperr <= idemtol) ) breakloop = 1;

#if REALSIZE == 4
    if ( sp2convint == 1 && fabs(limdiff) <= idemtol) breakloop = 1;
#elif REALSIZE == 8    
    if ( sp2convint == 1 && abs(limdiff) <= idemtol) breakloop = 1;
#endif  

  }
    
  if (iter == 100) {
    printf("SP2 pure has reached 100 iterations: something is wrong! \n");
    exit(1);
  }
  
  // x0 = 2.0 * x0

  M_MultiplyScalarMgpu( &TWO, x0);

  M_PullMgpu(x0);
  
  M_DeallocateLocal(xtmp);
  M_DeallocateDevice(xtmp);
  M_DeallocateDevice(x0);

}

void sp2pure_nospin4(REAL bndfil, int  hdim, REAL *x0_pointer, REAL maxeval,
             REAL *h_pointer, REAL maxminusmin, int minsp2iter,
             int sp2convint) {

  int iter, breakloop;
  REAL trx, trxtmp, trxold;
  REAL occ, limdiff;
  REAL idemperr = ZERO, idemperr1 = ZERO, idemperr2 = ZERO;
  REAL idemtol = SMALL_NUMBER;
  Matrix xtmp, x0;

  M_Init(xtmp, hdim, hdim);
  M_InitWithLocal(x0, x0_pointer, hdim, hdim);

  occ=bndfil*hdim;

  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      x0_pointer[i*hdim + j] = -h_pointer[i*hdim + j]/maxminusmin;
    }
  }

  for (int i = 0; i < hdim; i++) {
    x0_pointer[i*hdim + i] = maxeval/maxminusmin + x0_pointer[i*hdim + i];
  }

  // Copy x0 GPU
  M_Push( x0 );

  trx =  M_Trace( x0 );

  //fprintf(stdout, "initial trace x0 = %f \n", trx);

  iter = 0;

  breakloop = 0;

  while( breakloop == 0 && iter < 100 ) {

    iter++;

    M_Copy(x0, xtmp);

    // On multiple GPUs: xtmp = x0 * x0 + x0
    M_Multiply( &MINUS1, x0, x0, &ONE, xtmp);

    trxtmp = M_Trace( xtmp );

    trxold = trx;
         
#if REALSIZE == 4
    limdiff = fabs(trx - trxtmp - occ) - fabs(trx + trxtmp - occ);
#elif REALSIZE == 8
    limdiff = abs(trx - trxtmp - occ) - abs(trx + trxtmp - occ); 
#endif

    if (limdiff > idemtol ) {

      // x0 = x0 + xtmp
      M_MultiplyScalarSum( &ONE, xtmp, x0);
      
      trx = trx + trxtmp;

    } else if ( limdiff < -idemtol ) {
      
      // x0 = x0 - xtmp
      M_MultiplyScalarSum( &MINUS1, xtmp, x0);
      
      trx = trx - trxtmp;
      
    }
    
//    fprintf(stdout, "%d trace x0 = %f \n", iter, trx);

    idemperr2 = idemperr1;
    idemperr1 = idemperr;
#if REALSIZE == 4
    idemperr = fabs( trx - trxold );
    //idemperr = fabs( trxtmp );
#elif REALSIZE == 8
    idemperr = abs( trx - trxold );
    //idemperr = abs( trxtmp );
#endif

        //fprintf(stdout, "%d %g %g %g %g\n", iter, idemperr, trx, limdiff, trxtmp);

/*
    if ( sp2convint == 0 && iter >= minsp2iter && idemperr2 <= idemperr )
       breakloop = 1;
*/

    if ( sp2convint == 0 && iter >= minsp2iter &&
         (idemperr2 <= idemperr || idemperr <= idemtol) ) breakloop = 1;
    
#if REALSIZE == 4
    if ( sp2convint == 1 && fabs(limdiff) <= idemtol) breakloop = 1;
#elif REALSIZE == 8    
    if ( sp2convint == 1 && abs(limdiff) <= idemtol) breakloop = 1;
#endif

  }
    
  if (iter == 100) {
    printf("SP2 pure has reached 100 iterations: something is wrong! \n");
    exit(1);
  }

  // x0 = 2.0 * x0
  M_MultiplyScalar( &TWO, x0);
  M_Pull(x0);

  M_DeallocateLocal(xtmp); 
  M_DeallocateDevice(xtmp);
  M_DeallocateDevice(x0);
    
}

void sp2pure_spin3(REAL bndfil, int  hdim, REAL *rhoup_pointer, REAL *rhodown_pointer, REAL maxeval, REAL *hup_pointer, REAL *hdown_pointer, REAL maxminusmin, int minsp2iter, int sp2convint) {

  int iter, breakloop;
  REAL trx, totne, trxtmp, limdiff;
  REAL idemtol = SMALL_NUMBER;

  Matrix xtmpup, xtmpdown, rhoup, rhodown;

  M_Init(xtmpup, hdim, hdim);
  M_Init(xtmpdown, hdim, hdim);
  M_InitWithLocal(rhoup, rhoup_pointer, hdim, hdim);
  M_InitWithLocal(rhodown, rhodown_pointer, hdim, hdim);
  
  //
  // We're also using Niklasson's scheme to determine convergence
  //

  REAL idemperr = ZERO, idemperr1 = ZERO, idemperr2 = ZERO, trxold;

  totne = TWO*bndfil*hdim;

  //
  // Start by remapping the two H matrices such that 
  // all eigenvalues are [0:1]. We have the Gersgorin bounds
  // for both Hup and Hdown
  //

  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      rhoup_pointer[i*hdim + j] = -hup_pointer[i*hdim + j]/maxminusmin;
      rhodown_pointer[i*hdim + j] = -hdown_pointer[i*hdim + j]/maxminusmin;
    }
  }

  for (int i = 0; i < hdim; i++) {

    rhoup_pointer[i*hdim + i] = maxeval/maxminusmin + 
      rhoup_pointer[i*hdim + i];
    rhodown_pointer[i*hdim + i] = maxeval/maxminusmin + 
      rhodown_pointer[i*hdim + i];

  }

  M_Push(rhoup);
  M_Push(rhodown);

  trx = M_Trace(rhoup) + M_Trace(rhodown);

  iter = 0;
  breakloop = 0;

  while (breakloop == 0 && iter < 100) {

    iter++;

    M_Copy(rhoup, xtmpup);
    M_Copy(rhodown, xtmpdown);

    M_Multiply( &MINUS1, rhoup, rhoup, &ONE, xtmpup);
    M_Multiply( &MINUS1, rhodown, rhodown, &ONE, xtmpdown);
    
    trxtmp = M_Trace( xtmpup ) + M_Trace( xtmpdown );
    
#if REALSIZE == 4
    limdiff = fabs(trx - trxtmp - totne) - fabs(trx + trxtmp - totne);
#elif REALSIZE == 8
    limdiff = abs(trx - trxtmp - totne) - abs(trx + trxtmp - totne);
#endif

    trxold = trx;

    if (limdiff > idemtol ) {
      
      // rhoup = rhoup + xtmpup
      // rhodown = rhodown + xtmpdown
      M_MultiplyScalarSum( &ONE, xtmpup, rhoup);
      M_MultiplyScalarSum( &ONE, xtmpdown, rhodown);
      
      trx = trx + trxtmp;
      
    } else if ( limdiff < -idemtol ) {
      
      // rhoup = rhoup - xtmpup
      // rhodown = rhodown - xtmpdown
      M_MultiplyScalarSum( &MINUS1, xtmpup, rhoup);
      M_MultiplyScalarSum( &MINUS1, xtmpdown, rhodown);
      
      trx = trx - trxtmp;
      
    }
    
    idemperr2 = idemperr1;	
    idemperr1 = idemperr;
#if REALSIZE == 4
    idemperr = fabs( trx - trxold ); 
#elif REALSIZE == 8
    idemperr = abs( trx - trxold );
#endif

//    printf("%d %g %g\n", iter, idemperr, limdiff);
    
    if ( sp2convint == 0 && iter >= minsp2iter && 
	 (idemperr2 <= idemperr || idemperr <= idemtol) ) breakloop = 1;
    
#if REALSIZE == 4
    if ( sp2convint == 1 && fabs(limdiff) <= idemtol) breakloop = 1;
#elif REALSIZE == 8 
    if ( sp2convint == 1 && abs(limdiff) <= idemtol) breakloop = 1;
#endif

    cudaThreadSynchronize();

  }

  if (iter == 100) 
    printf("SP2 pure has reached 100 iterations: something is wrong! \n");
  
  M_Pull(rhoup);
  M_Pull(rhodown);

  M_DeallocateLocal(xtmpup);
  M_DeallocateDevice(xtmpup);
  M_DeallocateLocal(xtmpdown);
  M_DeallocateDevice(xtmpdown);
  M_DeallocateDevice(rhoup);
  M_DeallocateDevice(rhodown);

}
 
