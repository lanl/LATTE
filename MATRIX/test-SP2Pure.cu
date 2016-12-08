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

void sp2pure_nospin3(REAL bndfil, int  hdim, REAL *x0_pointer, REAL maxeval,
             REAL *h_pointer, REAL maxminusmin, int minsp2iter, 
	     int sp2convint) {

  int iter, breakloop;
  REAL trx, trxtmp, trxold;
  REAL occ, limdiff;
  REAL idemperr = 0.0, idemperr1 = 0.0, idemperr2 = 0.0;
  REAL idemtol = 1.0e-14;
  Matrix xtmp, x0, h;

  M_Init(xtmp, hdim, hdim);
  M_InitWithLocal(x0, x0_pointer, hdim, hdim);
  M_InitWithLocal(h, h_pointer, hdim, hdim);

  occ=bndfil*hdim;

  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      x0_pointer[i*hdim + j] = 0.0;
    }
  }

  trx = 0.0;

  for (int i = 0; i < hdim; i++) {
    x0_pointer[i*hdim + i] = 1.0;
    trx += x0_pointer[i*hdim + i];
  }

  /*  for (int i = 0; i < hdim; i++) {
    for (int j = 0; j < hdim; j++) {
      x0_pointer[i*hdim + j] = -h_pointer[i*hdim + j]/maxminusmin;
    }
  }

  trx = 0.0;

  for (int i = 0; i < hdim; i++) {
    x0_pointer[i*hdim + i] = maxeval/maxminusmin + x0_pointer[i*hdim + i];
    trx += x0_pointer[i*hdim + i];
  }

  */
  M_Push(x0);

  iter=0;

  breakloop = 0;

  printf("trx = %f \n", trx);

  M_Copy(x0, xtmp);
  cudaThreadSynchronize();
  trxtmp = M_Trace( xtmp );
  cudaThreadSynchronize();
  printf("trxtmp = %f \n", trxtmp);
  printf("hdim = %d \n", hdim);
  M_Add(x0, xtmp, x0);
  cudaThreadSynchronize();
  trx = M_Trace( x0 );
  cudaThreadSynchronize();
  printf("trx2 = %f \n", trx);

/*  while( breakloop == 0 && iter < 100 ) {
    
    iter++;
    
    M_Copy(x0, xtmp);

    M_Multiply(-1.0, x0, x0, 1.0, xtmp);

    trxtmp = M_Trace( xtmp );

    trxold = trx;
    
    limdiff = fabs(trx - trxtmp - occ) - fabs(trx + trxtmp - occ);
    
    if (limdiff > idemtol ) {
      
      M_Add(x0, xtmp, x0);
      
      trx = trx + trxtmp;
      
    } else if ( limdiff < -idemtol ) {
      
      M_Subtract(x0, xtmp, x0);
      
      trx = trx - trxtmp;
      
    }
    
    idemperr2 = idemperr1;	
    idemperr1 = idemperr;
    idemperr = fabs( trx - trxold );
    
    printf("%d %g %g %g %g\n", iter, idemperr, trx, limdiff, trxtmp);
    
    if ( sp2convint == 0 && iter >= minsp2iter && 
	 (idemperr2 <= idemperr || idemperr <= idemtol) ) breakloop = 1;
    
    if ( sp2convint == 1 && fabs(limdiff) <= idemtol) breakloop = 1;
    
  }
*/
  if (iter == 100) 
    printf("SP2 pure has reached 100 iterations: something is wrong! \n");

  
  M_Multiply(2.0, x0, x0);
  M_Pull(x0);
  
  M_DeallocateLocal(xtmp);
  M_DeallocateDevice(xtmp);
  M_DeallocateDevice(x0);
  M_DeallocateDevice(h);
}


void sp2pure_spin3(REAL bndfil, int  hdim, REAL *rhoup_pointer, REAL *rhodown_pointer, REAL maxeval, REAL *hup_pointer, REAL *hdown_pointer, REAL maxminusmin, int minsp2iter, int sp2convint) {

  int iter, breakloop;
  REAL trx, totne, trxtmp, limdiff;
  REAL idemtol = 1.0e-14;

  Matrix xtmpup, xtmpdown, rhoup, rhodown, hup, hdown;
  M_Init(xtmpup, hdim, hdim);
  M_Init(xtmpdown, hdim, hdim);
  M_InitWithLocal(rhoup, rhoup_pointer, hdim, hdim);
  M_InitWithLocal(rhodown, rhodown_pointer, hdim, hdim);
  M_InitWithLocal(hup, hup_pointer, hdim, hdim);
  M_InitWithLocal(hdown, hdown_pointer, hdim, hdim);
  
  //
  // We're also using Niklasson's scheme to determine convergence
  //

  REAL idemperr = 0.0, idemperr1 = 0.0, idemperr2 = 0.0, trxold;

  totne = 2.0*bndfil*hdim;


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

  trx = 0.0;

  for (int i = 0; i < hdim; i++) {

    rhoup_pointer[i*hdim + i] = maxeval/maxminusmin + 
      rhoup_pointer[i*hdim + i];
    rhodown_pointer[i*hdim + i] = maxeval/maxminusmin + 
      rhodown_pointer[i*hdim + i];

    trx += rhoup_pointer[i*hdim + i] + rhodown_pointer[i*hdim + i];
  }

  M_Push(rhoup);
  M_Push(rhodown);

  iter = 0;
  breakloop = 0;

  while (breakloop == 0 && iter < 100) {

    iter++;

    M_Copy(rhoup, xtmpup);
    M_Copy(rhodown, xtmpdown);

    M_Multiply(-1.0, rhoup, rhoup, 1.0, xtmpup);
    M_Multiply(-1.0, rhodown, rhodown, 1.0, xtmpdown);
    
    trxtmp = M_Trace( xtmpup ) + M_Trace( xtmpdown );
    
    limdiff = fabs(trx - trxtmp - totne) - fabs(trx + trxtmp - totne);

    trxold = trx;

    if (limdiff > idemtol ) {
      
      M_Add(rhoup, xtmpup, rhoup);
      M_Add(rhodown, xtmpdown, rhodown);
      
      trx = trx + trxtmp;
      
    } else if ( limdiff < -idemtol ) {
      
      M_Subtract(rhoup, xtmpup, rhoup);
      M_Subtract(rhodown, xtmpdown, rhodown);
      
      trx = trx - trxtmp;
      
    }
    
    idemperr2 = idemperr1;	
    idemperr1 = idemperr;
    idemperr = fabs( trx - trxold ); 

//    printf("%d %g %g\n", iter, idemperr, limdiff);
    
    if ( sp2convint == 0 && iter >= minsp2iter && 
	 (idemperr2 <= idemperr || idemperr <= idemtol) ) breakloop = 1;
    
    if ( sp2convint == 1 && fabs(limdiff) <= idemtol) breakloop = 1;

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
  M_DeallocateDevice(hup);
  M_DeallocateDevice(hdown);
}
 
