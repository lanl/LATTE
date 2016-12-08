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

void sp2fermi_nospin(REAL bndfil, int hdim, REAL *x0_pointer, 
                     REAL maxeval, REAL *h_pointer, REAL maxminusmin, 
                     REAL *chempot_pointer, int norecs, int *signlist_pointer,
                     REAL *beta0_pointer, REAL breaktol) {
  //
  // This is an implementation of Niklasson's SP2 algorithm for the
  // Fermi operator (i.e., finite temperature SP2)
  //
  // GERSHGORIN and SP2FERMIINIT must be run first to initialize 
  //

  int i, j, ii;
  int iter, breakloop;
  REAL occ = bndfil*hdim;
  REAL trx, trxomx;
  REAL lambda;
  REAL maxshift = ONE;
  REAL occerror;
  REAL preverror = ZERO, preverror2 = ZERO, preverror3 = ZERO;
  REAL gemm_alpha, gemm_beta;
  Matrix x0, xtmp;

  M_InitWithLocal(x0, x0_pointer, hdim, hdim);
  M_Init(xtmp, hdim, hdim); 

  iter = 0;
  breakloop = 0;

  int *signlist;

  signlist = (int*) malloc(norecs*sizeof(int));
  
  for (i = 0; i < norecs; i++) {
    signlist[ i ] = signlist_pointer[ i ];
  }     

  REAL chempot=*chempot_pointer;
  REAL beta0=*beta0_pointer;
  
  while (breakloop == 0 && iter < 100) {
    
    iter = iter + 1;
    
    for( i = 0; i < hdim; i++) {
      for(j = 0; j < hdim; j++) {
	
	x0_pointer[ i*hdim + j ] = -h_pointer[ i*hdim + j ]/maxminusmin;
	
      }
    }
    
    REAL buildfact;
    
    buildfact = (maxeval - chempot)/maxminusmin;
    
    for ( i = 0; i < hdim; i++ ) {
      
      x0_pointer[ i*hdim + i ] += buildfact ;
      
    }
    
    M_Push(x0);
    
    for(ii = 0; ii < norecs; ii++ ) {
      
      M_Copy( x0, xtmp );
      
      if (signlist[ii] > 0) {
	
	M_Multiply( &MINUS1, xtmp, xtmp, &TWO, x0) ;
	
      } else {
	
	M_Multiply( xtmp, xtmp, x0) ;
	
      }
      
    }

    // Trace ( X )
    
    trx = M_Trace( x0 );
    
    // Trace ( X^2 )

    REAL trx2;
    
    trx2 = M_TraceX2( x0 );

    // Trace ( X(1 - X) )

    trxomx = trx - trx2;

//	printf("%f %f \n", trx, trx2);

    lambda = (occ - trx)/(beta0 * trxomx);

    if (lambda > maxshift) lambda = maxshift;
    if (lambda < -maxshift) lambda = -maxshift;

    chempot = chempot + lambda;

    preverror3 = preverror2;
    preverror2 = preverror;
    preverror = occerror;
#if REALSIZE==4
    occerror = fabs(occ - trx);
#elif REALSIZE==8
    occerror = abs(occ - trx);
#endif

//    printf("%d %f %f \n ", iter, chempot, occerror);

    if (sizeof(REAL) == 8) {
      if (occerror < breaktol) {
	breakloop = 1;
      }
    }
    else {
      if (occerror == preverror || occerror == preverror2 || occerror == preverror3 || iter == 10 ) {
	breakloop = 1;
      }
    }
  }

  if (iter == 100) {
    printf("sp2fermi is not converging: stop!\n");
    exit(1);
  }

  gemm_alpha = MINUS2*lambda*beta0;
  gemm_beta = TWO - gemm_alpha;

  M_Copy( x0, xtmp );

  M_Multiply( &gemm_alpha, xtmp, xtmp, &gemm_beta, x0 );

  M_Pull(x0);

  *chempot_pointer=chempot;
  //  *kbt_pointer=kbt;

//  printf("%f \n", chempot);

  free(signlist);
  M_DeallocateDevice(x0);
  M_DeallocateDevice(xtmp);
  M_DeallocateLocal(xtmp);

}
/*
void sp2fermi_spin(REAL bndfil, int hdim, REAL *rhoup_ptr, REAL *rhodown_ptr, REAL maxeval, REAL *hup, REAL *hdown, REAL maxminusmin, REAL *chempot_pointer,
			  int norecs, REAL *kbt_pointer, REAL *beta0_pointer, REAL breaktol) {
  // This is an implementation of Niklasson's SP2 algorithm for the
  // Fermi operator (i.e., finite temperature SP2)
  // GERSHGORIN and SP2FERMIINIT must be run first to initialize 
  int i, j, ii;
  int iter=0, breakloop=0;
  REAL occ=2.0*bndfil*hdim, trx, trx1, lambda=0.0, maxshift = 1.0, occerror=0.0;
  REAL preverror = 0.0, preverror2 = 0.0, preverror3 = 0.0;
  REAL beta0=*beta0_pointer, chempot=*chempot_pointer;
  REAL totne=2.0*bndfil*hdim;
  
  Matrix rhoup, rhodown, x2up, x2down;
  M_InitWithLocal(rhoup, rhoup_ptr, hdim, hdim);
  M_InitWithLocal(rhodown, rhodown_ptr, hdim, hdim);
  M_Init(x2up, hdim, hdim);
  M_Init(x2down, hdim, hdim);
  
  while ( breakloop == 0 ) {
    iter++;
    if (iter == 50) {
      printf("sp2fermi is not converging: stop\n");
      exit(1);
    }
    
    for(i = 0; i<hdim; i++) {
      for (j = 0; j<hdim; j++) {
        if (i == j) {
          rhoup.Local[i+i*hdim] = (maxeval - hup[i+i*hdim] - chempot)/maxminusmin;
          rhodown.Local[i+i*hdim] = (maxeval - hdown[i+i*hdim] - chempot)/maxminusmin;
	}
        else {
          rhoup.Local[j+i*hdim] = -hup[j+i*hdim]/maxminusmin;
          rhoup.Local[i+j*hdim] = rhoup.Local[j+i*hdim];

          rhodown.Local[j+i*hdim] = -hdown[j+i*hdim]/maxminusmin;
          rhodown.Local[i+j*hdim] = rhodown.Local[j+i*hdim];
	}
      }
    }

    M_Push(rhoup);
    M_Push(rhodown);
    for (ii = 0; ii<norecs; ii++) {
      M_Multiply(rhoup, rhoup, x2up);
      M_Multiply(rhodown, rhodown, x2down);

      if (signlist[ii]==1) {
	M_Multiply(TWO, rhoup, rhoup);
	M_Subtract(rhoup, x2up, rhoup);
	M_Multiply(TWO, rhodown, rhodown);
	M_Subtract(rhodown, x2down, rhodown);
      }
      if (signlist[ii]==-1) {
	M_Copy(x2up, rhoup);
	M_Copy(x2down, rhodown);
      }
    }

    trx = M_Trace(rhoup) + M_Trace(rhodown);
    trx1 = beta0*(trx - (M_TraceX2(rhoup) + M_TraceX2(rhodown)));

    if (fabs(trx1)>1e-6)
      lambda = (totne - trx)/trx1;
    else {
      if (totne > trx)
        lambda = maxshift;
      else
        lambda = -maxshift;
    }

//printf(" lambda=%f\n", lambda);

    if (fabs(lambda) > maxshift)
      lambda = sign(maxshift, lambda);

    chempot = chempot + lambda;
//printf(" chemical potential=%f\n", chempot);

    preverror3 = preverror2;
    preverror2 = preverror;
    preverror = occerror;
        
    occerror = fabs(occ - trx);
//printf(" occupation error=%f\n", occerror);

    // how we figure if we've reached convergence. an absolute
    // tolerance works well in double precision, but in single
    // precision we need to look for noise when we're near 
    // self-consistency

    if (sizeof(REAL) == 8) {
      if (occerror < breaktol)
        breakloop = 1;
    }
    else {
      if (occerror == preverror || occerror == preverror2 || occerror == preverror3 || iter == 25 )
        breakloop = 1;
    }
  }
  M_Pull(rhoup);
  M_Pull(rhodown);

  M_DeallocateDevice(x2up);
  M_DeallocateLocal(x2up);
  M_DeallocateDevice(x2down);
  M_DeallocateLocal(x2down);
  M_DeallocateDevice(rhoup);
  M_DeallocateDevice(rhodown);

  *chempot_pointer=chempot;
}
*/
