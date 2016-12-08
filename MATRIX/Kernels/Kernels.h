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

#ifndef LATTE_KERNELS_H
#define LATTE_KERNELS_H

#if REALSIZE==4
  #define REAL float
#elif REALSIZE==8
  #define REAL double
#endif

#define GPU_ZERO (REAL)0.0
#define GPU_ONE (REAL)1.0

__global__ void MatrixAssembleKernel(int const size, int n, REAL *A, REAL *A2, int b, int sub);
__global__ void MatrixFastTraceX2Kernel(REAL *A, REAL *trace, int M, int num_threads);
__global__ void MatrixFastTraceX2Kernel(int const size, REAL *A, REAL *trace);
__global__ void FastDotProductKernel(int const size, REAL *A, REAL *B, REAL *dot);
__global__ void SubtractMatrixKernel(int const size, REAL *A, REAL *B, REAL *C);
__global__ void AddIdentityKernel(int const size, REAL *A);
__global__ void MultiplyScalarMatrixKernel(int const size, REAL Scalar, REAL *A, REAL *B);
__global__ void MultiplyScalarMatrixAddKernel(int const size, REAL Scalar, REAL *A, REAL Scalar2, REAL *B, REAL *C);
__global__ void MultiplyScalarMatrixSubKernel(int const size, REAL Scalar, REAL *A, REAL Scalar2, REAL *B, REAL *C);
__global__ void MultiplyScalarMatrixAddMatrixKernel(int const size, REAL Scalar, REAL *A, REAL *B, REAL *C);
__global__ void MultiplyScalarMatrixSubMatrixKernel(int const size, REAL Scalar, REAL *A, REAL *B, REAL *C);
__global__ void SubtractVectorKernel(int const size, REAL k, REAL *A, REAL *B, REAL *C);
__global__ void CGIterateKernel(int const M, REAL *p0, REAL *tmpmat, REAL *r0, REAL *bo, REAL *error2_ptr);
__global__ void AddMatrixKernel(int const size, REAL *A, REAL *B, REAL *C);
__global__ void MatrixFastTraceKernel(int const m, const int size, REAL *A, REAL *result, int offset);
__global__ void AddVectorKernel(int const size, REAL k, REAL *A, REAL *B, REAL *C);
__global__ void MatrixTraceKernel(REAL *A, REAL *trace, int m, int num_threads);
__global__ void SubCopyKernel(int const size, int const istart, REAL *A, REAL *B);
__global__ void SubCopyTransKernel(int const size, int const istart, int AM, int AN, REAL *A, REAL *B);
__global__ void TransKernel(int const size, int AM, int AN, REAL *A, REAL *B);
#endif
