!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Constructs and list of vectors in the reciprocal space.
!! \brief This routine constructs a list of reciprocal vectors needed
!! to parallelize the reciprocal Ewald summation. These vectors are used
!! as inputs for the COULOMBEWALD routine. Briefly, this will compute a set
!! of vectors \f$K = \{\vec{k} \in \Re^3:k = l\vec{b_1} + m\vec{b_2}
!! + n\vec{b_3},\, \mathrm{with}\, 0 \leq l \leq l_{max}, 0 \leq m \leq m_{max}
!! \, \mathrm{and}\, 0 \leq n \leq n_{max}\}\f$. In this case, \f$ \vec{b_1},
!! \vec{b_2},  \vec{b_3} \f$ are the reciprocal vectors previously
!! constructed in INITCOULOMB.
!!
!! \param RECIP_VECTORS Reciprocal lattice vectors.
!! \param K1_LIST (global output) list of Kx values.
!! \param K2_LIST (global output) list of Ky values.
!! \param K3_LIST (global output) list of Kz values.
!! \param KSQ_LIST list of K2 (k squared) values.
!! \param KCUTOFF (global input) Cutoff value for k.
!! \param LMAX (global input) Number of point in the b1 direction.
!! \param MMAX (global input) Number of point in the b2 direction.
!! \param NMAX (global input) Number of point in the b3 direction.
!! \param NK (global output) Number of total points in the reciprocal space.
!!
SUBROUTINE GET_K_LISTS(RECIP_VECTORS)

  USE CONSTANTS_MOD
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE
  INTEGER                            ::  I, L, M, N
  INTEGER                            ::  LMIN, MMIN, NMIN
  REAL(LATTEPREC)                    ::  L11, L12, L13, M21
  REAL(LATTEPREC)                    ::  M22, M23
  REAL(LATTEPREC)                    ::  K(3), K2
  REAL(LATTEPREC), INTENT(IN)        ::  RECIP_VECTORS(3,3)

  IF (EXISTERROR) RETURN

  !Maximum value of vectors in the reciprocal space
  NK = (2*LMAX+1)*(2*MMAX+1)*(2*NMAX+1)

  IF(.NOT.ALLOCATED(K1_LIST))ALLOCATE(K1_LIST(NK))
  IF(.NOT.ALLOCATED(K2_LIST))ALLOCATE(K2_LIST(NK))
  IF(.NOT.ALLOCATED(K3_LIST))ALLOCATE(K3_LIST(NK))
  IF(.NOT.ALLOCATED(KSQ_LIST))ALLOCATE(KSQ_LIST(NK))

  K1_LIST = 0.0
  K2_LIST = 0.0
  K3_LIST = 0.0
  KSQ_LIST = 0.0

  KCUTOFF2 = KCUTOFF * KCUTOFF

  I = 1
  LMIN = 0

  DO L = LMIN, LMAX

     IF (L == 0) THEN
        MMIN = 0
     ELSE
        MMIN = -MMAX
     ENDIF

     L11 = REAL(L,LATTEPREC)*RECIP_VECTORS(1,1)
     L12 = REAL(L,LATTEPREC)*RECIP_VECTORS(1,2)
     L13 = REAL(L,LATTEPREC)*RECIP_VECTORS(1,3)

     DO M = MMIN, MMAX

        NMIN = -NMAX

        IF(L .EQ. 0 .AND. M .EQ. 0) NMIN = 1

        M21 = L11 + REAL(M,LATTEPREC)*RECIP_VECTORS(2,1)
        M22 = L12 + REAL(M,LATTEPREC)*RECIP_VECTORS(2,2)
        M23 = L13 + REAL(M,LATTEPREC)*RECIP_VECTORS(2,3)

        DO N = NMIN, NMAX

           K(1) = M21 + REAL(N,LATTEPREC)*RECIP_VECTORS(3,1)
           K(2) = M22 + REAL(N,LATTEPREC)*RECIP_VECTORS(3,2)
           K(3) = M23 + REAL(N,LATTEPREC)*RECIP_VECTORS(3,3)

           K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)

           IF (K2 .LE. KCUTOFF2) THEN

              K1_LIST(I) = K(1)
              K2_LIST(I) = K(2)
              K3_LIST(I) = K(3)
              KSQ_LIST(I) = K2
              I = I + 1

           END IF

        END DO

     END DO

  END DO

  NK = I-1

END SUBROUTINE GET_K_LISTS
