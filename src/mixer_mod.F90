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

!> To implement mixing schemes from the progress library
!!
MODULE MIXER_MOD


  USE CONSTANTS_MOD
  USE MYPRECISION
  USE COULOMBARRAY
  USE SETUPARRAY
  USE NONOARRAY
  USE DIAGARRAY    ! CHANGE ANDERS_CHANGE
  USE XBOARRAY     ! CHANGE ANDERS_CHANGE
  USE DMARRAY     ! CHANGE ANDERS
#ifdef PROGRESSON
  USE LATTEPARSER
  USE PRG_PULAYMIXER_MOD
  USE BML
  USE NONOARRAYPROGRESS
#endif
  PRIVATE

#ifdef PROGRESSON
  PUBLIC :: QMIXPRG
#endif

  !PUBLIC :: KERNELMIXER  ! this mixer is not used 
  PUBLIC :: KERNELPROPAGATION, FULLKERNELPROPAGATION
  PUBLIC :: PRECONDKERNELPROPAGATION, ADAPTPRECONDKERNEL
  PUBLIC :: DMKERNELMIXER, DMKERNELPROPAGATION, dP2MD, dP2MIXER

  !For mixing scheme
  LOGICAL, PUBLIC                      ::  MIXINIT = .FALSE.
  REAL(LATTEPREC), ALLOCATABLE, PUBLIC ::  DQIN(:,:), DQOUT(:,:)
  REAL(LATTEPREC), ALLOCATABLE, PUBLIC ::  QTMP1(:), QTMP2(:)
  REAL(LATTEPREC), PUBLIC              ::  SCFERROR
#ifdef PROGRESSON
  TYPE(MX_TYPE), PUBLIC                ::  MX
#endif

CONTAINS

  SUBROUTINE FULLKERNELPROPAGATION(MDITER)
    INTEGER, INTENT(IN) :: MDITER
    INTEGER :: I, J, N, ii, jj

    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:)
    !
    ALLOCATE(Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS))
    ALLOCATE(du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM))

    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
    write(6,*) 'MDITER', MDITER, norm2(Res) / SQRT(DBLE(NATS)) 

!    if (norm2(Res) <= 0.000001D0) then
!    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
!     dn2dt2 = MDMIX*Res
!    else
!    if (MDITER <= 20) then
    !    MDMIX = 1.D0
        dr = 0.D0*Res
        do I = 1,NATS !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
           dr(I) = 1.D0
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
              vi(:,I) = vi(:,I)/norm2(vi(:,I))
           enddo
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv
           ! ri(:,I) = dr + ((1.D0 - MDMIX)/MDMIX)*vi(:,I)
           ri(:,I) = dr 
           !du(:) = -MDMIX*ri(:,I)
           !wi(:,I) = -MDMIX*vi(:,I)
           du(:) = -ri(:,I)
           wi(:,I) = -vi(:,I)
           do J = 1,I-1
              du(:) = du(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
              wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
           enddo
           ui(:,I) = du/(1.D0 + dot_product(vi(:,I),du))
           dr(I) = 0.D0
        enddo
        FULL_K = 0.D0*FULL_K
        do I = 1,NATS 
           !FULL_K(I,I) = MDMIX
           FULL_K(I,I) = 1.0D0
        enddo
        call DGEMM('N','T',NATS,NATS,NATS,1.D0,ui,NATS,wi,NATS,1.D0,FULL_K,NATS)
!    endif
        ! update q corresponding to q = q - MATMUL(KK,Res)
        !DELTAQ = OLDDELTAQS + QMIX*Res
        dn2dt2 = MATMUL(FULL_K,Res)
!!        write(*,*) ' dn2dt2 FULL_K = ', dn2dt2(1:3)
!!        dn2dt2 = MDMIX*Res
!!        do I = 1,NATS  !! Let the approximate kernel act on the residual by individual rank-1 updates
!!           !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
!!           dn2dt2 = dn2dt2 + dot_product(wi(:,I),Res)*ui(:,I)
!!        enddo
!!        write(*,*) ' dn2dt2 Rank-Nats = ', dn2dt2(1:3)
!!        write(*,*) ' ------------------ '
!    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE)

  END SUBROUTINE FULLKERNELPROPAGATION

  SUBROUTINE PRECONDKERNELPROPAGATION(MDITER,LL)
    INTEGER, INTENT(IN) :: MDITER,LL
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Nocc, beta, eps
    INTEGER         :: INFO

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FO(:,:), FM(:,:), RI_T(:,:), WORK(:)
    INTEGER,         ALLOCATABLE :: IPIV(:)
    !
    ALLOCATE(Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS))
    ALLOCATE(du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(FO(LL,LL), FM(LL,LL), ri_t(LL,NATS),WORK(LL+LL*LL),IPIV(LL))

    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
    write(6,*) 'MDITER', MDITER, norm2(Res) / SQRT(DBLE(NATS)) 

!    if (norm2(Res) <= 0.00001D0) then
!    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
!     dn2dt2 = MDMIX*Res
!    else
        Res = MATMUL(FULL_K,Res) !! FULL_KK is the preconditioner
        dr = Res
        do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means preconditioning only!
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
           enddo
           vi(:,I) = vi(:,I)/norm2(vi(:,I))
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ ! v_{i+1} = paritla f(n+\lambda v_i} \partial lambda, eq 42
           dr = dq_dv - v
           dr = MATMUL(FULL_K,dr)
           ri(:,I) = dr
        enddo
        ri_t = transpose(ri)
        FO = MATMUL(ri_t,ri)
        FM = FO
        call DGETRF(LL, LL, FM, LL, IPIV, INFO)
        call DGETRI(LL, FM, LL, IPIV, WORK, LL+LL*LL, INFO)
        FO = MATMUL(FM,FO)
        dn2dt2 = 0.D0*Res
        do I = 1,LL
        do J = 1,LL
          dn2dt2 = dn2dt2 - FM(I,J)*dot_product(ri(:,J),Res)*vi(:,I)
        enddo
        enddo
!        dn2dt2 = -Res
!        do I = 1,LL
!        do J = 1,LL
!          dn2dt2 = dn2dt2 - FM(I,J)*dot_product(ri(:,J),Res)*(vi(:,I)-ri(:,J))
!        enddo
!        enddo
!    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,FO, FM, ri_t,WORK,IPIV)
  END SUBROUTINE PRECONDKERNELPROPAGATION

  SUBROUTINE ADAPTPRECONDKERNEL(MDITER,LL)
    implicit none
    INTEGER, INTENT(IN) :: MDITER,LL
    INTEGER :: I, J, N, II, JJ, K

    !REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
    !REAL(LATTEPREC) :: dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
    !REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS)
    !REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: Nocc, beta, eps
    !REAL(LATTEPREC) :: ORTHOH_SAVE(HDIM,HDIM)
    !REAL(LATTEPREC) :: ri_t(LL,NATS), IDENTRES(NATS)

    REAL(LATTEPREC) :: RESNORM, FEL, DTMP
    INTEGER         :: INFO, RANK
    REAL(LATTEPREC), ALLOCATABLE :: FO(:,:), FM(:,:),WORK(:)
    INTEGER,         ALLOCATABLE :: IPIV(:)

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: RI_T(:,:), IDENTRES(:)
    !
    ALLOCATE(Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS))
    ALLOCATE(du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(ri_t(LL,NATS),WORK(LL+LL*LL),IDENTRES(NATS))
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)

    RESNORM = NORM2(RES) / SQRT(DBLE(NATS)) 
   
    write(6,*) 'MDITER', MDITER, RESNORM
    IF (MDITER == 1) THEN
       CALL FULLKERNELPROPAGATION(MDITER)
    ELSE
        Res = MATMUL(FULL_K,Res) !! FULL_KK is the preconditioner
        dr = Res

        I = 0
        FEL = 1.D0
        DO WHILE (FEL > KERNELTOL)  !! LL is the number of rank-1 updates  LL = 0 means preconditioning only!
           I = I + 1
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
           enddo
           vi(:,I) = vi(:,I)/norm2(vi(:,I))
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv - v
           dr = MATMUL(FULL_K,dr)
           ri(:,I) = dr

           RANK = I
           ALLOCATE(FO(RANK, RANK), FM(RANK, RANK), IPIV(RANK))

           DO J = 1, RANK
              DO K = 1, RANK
                FO(K,J) = DOT_PRODUCT(RI(:,K), RI(:,J))
              ENDDO
           ENDDO

           CALL DGETRF(RANK, RANK, FO, RANK, IPIV, INFO)
           CALL DGETRI(RANK, FO, RANK, IPIV, WORK, LL+LL*LL, INFO)

           FM = FO
           DN2DT2 = 0.D0
           IDENTRES = 0.D0
           DO K = 1,RANK
           DO J = 1,RANK
             DTMP = FM(K,J)*dot_product(RI(:,J),RES)
             IdentRes = IdentRes + DTMP*RI(:,K)
             dn2dt2 = dn2dt2 - DTMP*VI(:,K)
           ENDDO
           ENDDO
           FEL = NORM2(IDENTRES - RES) / NORM2(RES)
           DEALLOCATE(FO, FM, IPIV)
        ENDDO 

        IF (RANK == LL ) THEN
           CALL FULLKERNELPROPAGATION(MDITER)
        ENDIF
    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,ri_t,WORK,IDENTRES)
  END SUBROUTINE ADAPTPRECONDKERNEL


  SUBROUTINE KERNELPROPAGATION(MDITER,LL)
    INTEGER, INTENT(IN) :: MDITER,LL
    INTEGER :: I, J, N, ii, jj

    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:)

    ALLOCATE(Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS))
    ALLOCATE(du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM))
    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0

    Res = DELTAQ - PNK(1,:)
!   if (norm2(Res) >= 0.00001D0) then
    if  (MDITER <= 0) then  !! typical choice <= 1, for really really hard cases <= 20
     dn2dt2 = MDMIX*Res
    else
        dr = Res
        do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means linear mixing
           vi(:,I) = dr/norm2(dr)
           do J = 1,I-1
              vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
              vi(:,I) = vi(:,I)/norm2(vi(:,I))
           enddo
           v(:) = vi(:,I)
           !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
           dq_dv = ZERO
           dq_v = v/norm2(v)
           DELTAQ = dq_v
           call coulombrspace
           call coulombewald
           call addqdep
           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
           BO = 2.D0*BO
           call deorthomyrho
           call getdeltaq_resp
           dq_dv = DELTAQ
           dr = dq_dv
           ri(:,I) = dr + ((1.D0 - MDMIX)/MDMIX)*vi(:,I)
           du(:) = -MDMIX*ri(:,I)
           wi(:,I) = -MDMIX*vi(:,I)
           do J = 1,I-1
              du(:) = du(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
              wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
           enddo
           ui(:,I) = du/(1.D0 + dot_product(vi(:,I),du))
        enddo
        ! update q corresponding to q = q - MATMUL(KK,Res)
        !DELTAQ = OLDDELTAQS + QMIX*Res
        dn2dt2 = MDMIX*Res
        do I = 1,LL  !! Let the approximate kernel act on the residual by individual rank-1 updates
           !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
           dn2dt2 = dn2dt2 + dot_product(wi(:,I),Res)*ui(:,I)
        enddo
    endif
    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE
    DELTAQ = DELTAQ_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE)
  END SUBROUTINE KERNELPROPAGATION

  !SUBROUTINE KERNELMIXER(PITER,LL)
  !  INTEGER, INTENT(IN) :: PITER,LL
  !  INTEGER :: I, J, N, ii, jj
  !  REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
  !  REAL(LATTEPREC) :: du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
  !  REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS)
  !  REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
  !  REAL(LATTEPREC) :: Nocc, beta, eps, FCOUL_SAVE(3,NATS)
  !  !REAL(LATTEPREC) :: ev(HDIM), mu0 , fe(HDIM), DDT(HDIM,HDIM), mu1, trX, trDDT
  !  REAL(LATTEPREC) :: ORTHOH_SAVE(HDIM,HDIM)
  !  !
  !  DELTAQ_SAVE = DELTAQ
  !  COULOMBV_SAVE = COULOMBV
  !  FCOUL_SAVE = FCOUL
  !  BO_SAVE = BO
  !  ORTHOH_SAVE = ORTHOH
  !  H_SAVE = H
  !  H_0 = H0
  !  H0 = 0.D0

  !  Res = DELTAQ - OLDDELTAQS
  !  !if  (PITER <= 10) then  !! typical choice <= 1, for really really hard cases <= 20
  !  if ((norm2(Res)/sqrt(1.D0*NATS) > 1e-1).OR.(PITER<=1)) then
  !    DELTAQ = OLDDELTAQS + QMIX*Res
  !  else
  !      dr = Res
  !      do I = 1,LL !! LL is the number of rank-1 updates  LL = 0 means linear mixing
  !         vi(:,I) = dr/norm2(dr)
  !         do J = 1,I-1
  !            vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
  !            vi(:,I) = vi(:,I)/norm2(vi(:,I))
  !         enddo
  !         v(:) = vi(:,I)
  !         !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
  !         dq_dv = ZERO
  !         dq_v = v/norm2(v)
  !         DELTAQ = dq_v
  !         call coulombrspace
  !         call coulombewald
  !         call addqdep
  !         call orthomyh
  !         Nocc = BNDFIL*float(HDIM)
  !         beta = 1.D0/KBT
  !         !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
  !         call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
  !         BO = 2.D0*BO
  !         call deorthomyrho
  !         call getdeltaq_resp
  !         dq_dv = DELTAQ
  !         dr = dq_dv
  !         ri(:,I) = dr + ((1.D0 - QMIX)/QMIX)*vi(:,I)
  !         du(:) = -QMIX*ri(:,I)
  !         wi(:,I) = -QMIX*vi(:,I)
  !         do J = 1,I-1
  !            du(:) = du(:) - dot_product(wi(:,J),ri(:,I))*ui(:,J)
  !            wi(:,I) = wi(:,I) - dot_product(ui(:,J),vi(:,I))*wi(:,J)
  !         enddo
  !         ui(:,I) = du/(1.D0 + dot_product(vi(:,I),du))
  !      enddo
  !      ! update q corresponding to q = q - MATMUL(KK,Res)
  !      DELTAQ = OLDDELTAQS + QMIX*Res
  !      do I = 1,LL  !! Let the approximate kernel act on the residual by individual rank-1 updates
  !         DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
  !      enddo
  !  endif
  !  COULOMBV = COULOMBV_SAVE
  !  BO = BO_SAVE
  !  H0 = H_0
  !  H = H_SAVE
  !  FCOUL = FCOUL_SAVE
  !  ORTHOH = ORTHOH_SAVE
  !END SUBROUTINE KERNELMIXER        

!!! ANDERS CHANGE ADD NEW SUBROUTINE
  SUBROUTINE DMKERNELPROPAGATION(PITER)
    INTEGER, INTENT(IN) :: PITER
    INTEGER :: I, J, N, ii, jj
    REAL(LATTEPREC) :: Nocc, beta, eps,ndDO_dU,nDO

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), FCOUL_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: X(:,:), YY(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:), nDelta_DO(:,:), Delta_DS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_1(:,:),H1(:,:), dU(:,:)

    ALLOCATE(Res(NATS), dr(NATS), FCOUL_SAVE(3,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM), nDelta_DO(HDIM,HDIM),Delta_DS(HDIM,HDIM))
    ALLOCATE(H_1(HDIM,HDIM),H1(HDIM,HDIM), dU(HDIM,HDIM))
!
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0  ! SUCH THAT ADDQDEP ONLY INCLUDES RESPONSE PART
!
    Delta_DO = DOrth-PO_0
    nDO = 0.D0
    do I = 1, HDIM
      nDO = nDO + dot_product(Delta_DO(:,I),Delta_DO(:,I))
    enddo
    nDO =sqrt(nDO)  ! Twice the nDO in Developers version
    nDelta_DO = Delta_DO/nDO
!
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,Delta_DO,HDIM,ZERO,YY,HDIM)
    call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)
!
    Res = 0.D0
    do I = 1, NATS
      do K = H_INDEX_START(I), H_INDEX_END(I)
        Res(I) = Res(I) + Delta_DS(K,K)
      enddo
    enddo
!
    do I = 1, HDIM
    do J = 1, HDIM
       SU(I,J) = SMAT(I,J)*DFTB_U(J)
    enddo
    enddo
!
    call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,Delta_DS,HDIM,SU,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,Delta_DS,HDIM,ZERO,YY,HDIM)
    H1 = -0.25D0*(X+YY)/nDO   
    do I = 1,HDIM
    do J = 1,HDIM
       H_1(I,J) = H1(I,J) + H1(J,I)  ! Response in H from the Hubbard energy term
    enddo
    enddo
    DELTAQ = 2.D0*Res/nDO
    call coulombrspace
    call coulombewald
    call addqdep
    H = H + H_1  ! Total perturbation
    call orthomyh   ! ORTHOH is now the perturbation
    Nocc = BNDFIL*float(HDIM)
    beta = 1.D0/KBT
    !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
    call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
    !!! BO is now the DM response matrix with respect to the residual perturbation

    dU = BO + ((1-QMIX)/QMIX)*nDelta_DO  

    ndDO_dU = 0.D0
    do I = 1, HDIM
      ndDO_dU = ndDO_dU + dot_product(nDelta_DO(:,I),dU(:,I))
    enddo

    d2PO = QMIX*Delta_DO + (QMIX*QMIX*nDO/(1.D0-QMIX*ndDO_dU))*dU

    COULOMBV = COULOMBV_SAVE
    BO = BO_SAVE  
    DELTAQ = DELTAQ_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE

    DEALLOCATE(RES, DR, FCOUL_SAVE, DELTAQ_SAVE, COULOMBV_SAVE)
    DEALLOCATE(BO_SAVE, DELTA_DS, X, YY, ORTHOH_SAVE, DELTA_DO, NDELTA_DO)
    DEALLOCATE(H_0, H_1, H1, DU, H_SAVE)
  END SUBROUTINE DMKERNELPROPAGATION

!!! ANDERS CHANGE END NEW SUBROUTINE


  SUBROUTINE DMKERNELMIXER(PITER,MAXDQ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    REAL(LATTEPREC), INTENT(IN) :: MAXDQ
    !
    INTEGER :: I, J, N, ii, jj
    INTEGER :: K
    REAL(LATTEPREC) :: Nocc, beta, eps, ndDO_DU, nDO
    !
    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), FCOUL_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: X(:,:), YY(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:), nDelta_DO(:,:), Delta_DS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_1(:,:),H1(:,:), dU(:,:)

    ALLOCATE(Res(NATS), dr(NATS), FCOUL_SAVE(3,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM), nDelta_DO(HDIM,HDIM),Delta_DS(HDIM,HDIM))
    ALLOCATE(H_1(HDIM,HDIM),H1(HDIM,HDIM), dU(HDIM,HDIM))
!    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
!    BO_SAVE = BO
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    H0 = 0.D0  ! SUCH THAT ADDQDEP ONLY INCLUDES RESPONSE PART
!
    Delta_DO = DOrth-DOrth_old
    !write(*,*) ' Delta_DO = ', Delta_DO(1,:)/2.D0
    nDO = 0.D0
    do I = 1, HDIM
      nDO = nDO + dot_product(Delta_DO(:,I),Delta_DO(:,I))
    enddo
    nDO =sqrt(nDO)
    nDelta_DO = Delta_DO/nDO
!
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,Delta_DO,HDIM,ZERO,YY,HDIM)
    call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)
!
    do I = 1, NATS
      Res(I) = 0.D0
      do K = H_INDEX_START(I), H_INDEX_END(I)
        Res(I) = Res(I) + Delta_DS(K,K)
      enddo
    enddo
    !write(*,*) ' Res = ', Res(:)
!
    do I = 1, HDIM
    do J = 1, HDIM
       SU(I,J) = SMAT(I,J)*DFTB_U(J)
    enddo
    enddo
!
    call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,Delta_DS,HDIM,SU,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,Delta_DS,HDIM,ZERO,YY,HDIM)
    H1 = -0.25D0*(X+YY)/nDO   !! 
    do I = 1,HDIM
    do J = 1,HDIM
       H_1(I,J) = H1(I,J) + H1(J,I)  ! Response in H from the Hubbard energy term
    enddo
    enddo
    !write(*,*) ' H_1 = ', H_1(1,:)
    !write(*,*) ' nDO = ', nDO

    DELTAQ = 2.D0*Res/nDO
    call coulombrspace
    call coulombewald
    call addqdep
    H = H + H_1  ! Total perturbation
    call orthomyh
    Nocc = BNDFIL*float(HDIM)
    beta = 1.D0/KBT
    !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
    call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)

    dU = BO + ((1-QMIX)/QMIX)*nDelta_DO  ! Maybe 2*BO?

    ndDO_dU = 0.D0
    do I = 1, HDIM
      ndDO_dU = ndDO_dU + dot_product(nDelta_DO(:,I),dU(:,I))
    enddo

!    if (MAXDQ > 0.1D0) then
    if (PITER <= 10) then
      write(*,*) ' LINEAR MIX MIXER_MOD'
      BO = DOrth_old + QMIX*Delta_DO
    else
      write(*,*) ' RESPONSE MIX MIXER_MOD'
      BO = DOrth_old + QMIX*Delta_DO + (QMIX*QMIX*nDO/(1.D0-QMIX*ndDO_dU))*dU
    endif

    COULOMBV = COULOMBV_SAVE
!    BO = BO_SAVE  ! NOTE BO Has been replaced by the dD/dDelta_DO
    DELTAQ = DELTAQ_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE

    DEALLOCATE(RES, DR, FCOUL_SAVE, DELTAQ_SAVE, COULOMBV_SAVE)
    DEALLOCATE(BO_SAVE, DELTA_DS, X, YY, ORTHOH_SAVE, DELTA_DO, NDELTA_DO)
    DEALLOCATE(H_0, H_1, H1, DU, H_SAVE)
  END SUBROUTINE DMKERNELMIXER
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE dP2MD(PITER)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    !
    INTEGER :: I, J, K, N, ii, jj
    REAL(LATTEPREC), ALLOCATABLE :: Res(:), BO_SAVE(:,:), X(:,:), YY(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:), nDelta_DO(:,:),Delta_DS(:,:)

    ALLOCATE(Res(NATS),BO_SAVE(HDIM,HDIM),X(HDIM,HDIM), YY(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM), nDelta_DO(HDIM,HDIM),Delta_DS(HDIM,HDIM))
    
!
    BO_SAVE = BO

    Delta_DO = DOrth-PO_0

    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,Delta_DO,HDIM,ZERO,YY,HDIM)
    call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)

    do I = 1, NATS
      Res(I) = 0.D0
      do K = H_INDEX_START(I), H_INDEX_END(I)
        Res(I) = Res(I) + Delta_DS(K,K)
      enddo
    enddo
    !write(*,*) ' Res = ', Res(:)
    write(6,*) 'MDITER', PITER,NORM2(RES) / SQRT(DBLE(NATS)) 

    !write(6,*) 'maxval(dDO/dP) called from dP2MD'
    IF (PITER <= 1) then
      write(*,*) ' LINEAR MIX MIXER_MOD (MD)'
      d2PO = Delta_DO
      !write(6,*) PO_0
    ELSE
      write(*,*) ' RESPONSE MIX MIXER_MOD (MD)'
      call dP2(Delta_DO, d2PO)
    ENDIF

    BO = BO_SAVE

    DEALLOCATE(RES, BO_SAVE, X, YY, DELTA_DO, NDELTA_DO, DELTA_DS)
  END SUBROUTINE dP2MD
!
  SUBROUTINE dP2MIXER(PITER,SCF_ERR,MAXDQ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    REAL(LATTEPREC), INTENT(IN) :: SCF_ERR
    REAL(LATTEPREC), INTENT(IN) :: MAXDQ
    !
    INTEGER :: I, J, K, N, ii, jj

    REAL(LATTEPREC), ALLOCATABLE :: DU(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:)

    ALLOCATE(DU(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM))
!
    Delta_DO = DOrth-DOrth_old
    !nDO = 0.D0
    !do I = 1, HDIM
    !  nDO = nDO + dot_product(Delta_DO(:,I),Delta_DO(:,I))
    !enddo
    !nDO =sqrt(nDO)
    !nDelta_DO = Delta_DO/nDO
!
    !call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,Delta_DO,HDIM,ZERO,YY,HDIM)
    !call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
    !call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)

    !do I = 1, NATS
    !  Res(I) = 0.D0
    !  do K = H_INDEX_START(I), H_INDEX_END(I)
    !    Res(I) = Res(I) + Delta_DS(K,K)
    !  enddo
    !enddo
    !write(*,*) ' Res = ', Res(:)

    !IF(SCF_ERR > 0.01D0) THEN
    IF(SCF_ERR > 0.1D0) THEN
      !write(*,*) ' LINEAR MIX MIXER_MOD', QMIX, SCF_ERR
      dU = QMIX*Delta_DO
    ELSE
      !write(*,*) ' RESPONSE MIX MIXER_MOD', SCF_ERR
      call dP2(Delta_DO, dU) 
      !dU = QMIX*Delta_DO
    ENDIF
    BO = DOrth_old + dU

    DEALLOCATE(DU,DELTA_DO)

  END SUBROUTINE dP2MIXER



  SUBROUTINE dP2(Delta_DO,dP)
    IMPLICIT NONE
    REAL(LATTEPREC), intent(in)  :: Delta_DO(HDIM,HDIM)
    REAL(LATTEPREC), intent(out) :: dP(HDIM,HDIM)

    !
    INTEGER :: I, J, N, ii, jj
    INTEGER :: KI, MI, K, L
    REAL(LATTEPREC) :: ErrDM, TMP, Nocc, beta, eps, ErrDMn

    REAL(LATTEPREC), ALLOCATABLE :: IdentResDM(:,:), Res(:), FCOUL_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DM_V(:,:,:), UU(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: X(:,:), YY(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: nDelta_DO(:,:), Delta_DS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_1(:,:),H1(:,:), dU(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: MMX(:,:), MMX_I(:,:)

    ALLOCATE(IdentResDM(HDIM,HDIM), Res(NATS), FCOUL_SAVE(3,NATS))
    ALLOCATE(DM_V(HDIM,HDIM,NRANK), UU(HDIM,HDIM,NRANK))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(X(HDIM,HDIM), YY(HDIM,HDIM), ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(nDelta_DO(HDIM,HDIM), Delta_DS(HDIM,HDIM))
    ALLOCATE(H_1(HDIM,HDIM),H1(HDIM,HDIM), dU(HDIM,HDIM))
!    !
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    ORTHOH_SAVE = ORTHOH
    H_SAVE = H
    H_0 = H0
    ErrDM = 100_LATTEPREC
    eps   = 1.d-7

    dU = Delta_DO
    ! adpative rank-m, no preconditioning

    do I = 1, HDIM
    do J = 1, HDIM
       SU(I,J) = SMAT(I,J)*DFTB_U(J)
    enddo
    enddo

    K = 0
    ALLOCATE(MMX(1,1),MMX_I(1,1))
    do while ( ErrDM > KERNELTOL .AND. K < NRANK) 
       DEALLOCATE(MMX,MMX_I)
       H0 = 0.D0  ! SUCH THAT ADDQDEP ONLY INCLUDES RESPONSE PART
       K = K + 1
       DM_V(:,:,K) = dU(:,:)

       DO J = 1, K-1
          TMP = 0.0_LATTEPREC
          DO I = 1, HDIM
          DO L = 1, HDIM
            TMP = TMP + DM_V(L,I,K) * DM_V(L,I,J)
          ENDDO
          ENDDO
          DM_V(:,:,K) = DM_V(:,:,K) - TMP * DM_V(:,:,J)
       ENDDO
       
       TMP = NORM2(DM_V(:,:,K))
       IF (TMP<1.d-20) TMP = TMP + 1.d-20

      
       IF (TMP < 1.D-12) write(6,*) 'TMP is ZERO in dP2!'
       DM_V(:,:,K) = DM_V(:,:,K)/TMP
       !DM_V(:,:,K) = DM_V(:,:,K)/sqrt(TMP)

       call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,DM_V(:,:,K),HDIM,ZERO,YY,HDIM)
       call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
       call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)

       call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,Delta_DS,HDIM,SU,HDIM,ZERO,X,HDIM)
       call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,Delta_DS,HDIM,ZERO,YY,HDIM)

       H1 = -0.25D0*(X+YY)
       do I = 1,HDIM
       do J = 1,HDIM
          H_1(I,J) = H1(I,J) + H1(J,I)  ! Response in H from the Hubbard energy term
       enddo
       enddo

       DELTAQ = 0.0
       do I = 1,NATS
         do J = H_INDEX_START(I),H_INDEX_END(I)
          DELTAQ(I) = DELTAQ(I) + 2*Delta_DS(J,J)
        enddo
       enddo
   
       call coulombrspace
       call coulombewald
       call addqdep
       H = H + H_1  ! Total perturbation
       call orthomyh
       Nocc = BNDFIL*float(HDIM)
       beta = 1.D0/KBT

       ! call canonical response
       call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
       !call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)

       UU(:,:,K) = BO(:,:) - DM_V(:,:,K)
       dU(:,:) = UU(:,:,K)

       !---------------
       !check acccuracy
       !---------------
       ALLOCATE(MMX(K,K),MMX_I(K,K))
       do I = 1,K
       do J = 1,K
         TMP = 0.D0
         do L = 1,HDIM
         do KI = 1,HDIM
            TMP = TMP + UU(KI,L,I)*UU(KI,L,J)
         enddo
         enddo
         MMX(I,J) = TMP
       enddo
       enddo

       call Invert(MMX,MMX_I,K)
 
       IdentResDM = 0.D0
       do I = 1,K
       do J = 1,K
         TMP = 0.D0
         do L = 1,HDIM
         do KI = 1,HDIM
            TMP = TMP + UU(KI,L,J)*Delta_DO(KI,L)
         enddo
         enddo
         IdentResDM(:,:) = IdentResDM(:,:) + MMX_I(I,J)*TMP*UU(:,:,I)
       enddo
       enddo
       ErrDM = 0.D0
       ErrDMn = 0.D0
       do I = 1, HDIM
       do J = 1, HDIM
         ErrDM = ErrDM + (IdentResDM(J,I) - Delta_DO(J,I))**2
         ErrDMn = ErrDMn + (IdentResDM(J,I))**2
       enddo
       enddo
       IF(ErrDMn<1.d-20) ErrDMn = ErrDMn + 1.d-20
 
       ErrDM = sqrt(ErrDM)/sqrt(ErrDMn) ! Relative Error
       !write(*,*) ' ErrDM = ', ErrDM  !! Should go down as a function of K
       !-----
    enddo

    dP = 0.D0
    do I = 1,K
    do J = 1,K
       TMP = 0.D0
       do L = 1,HDIM
       do MI = 1,HDIM
          TMP = TMP + UU(MI,L,J)*Delta_DO(MI,L)
       enddo
       enddo
       dP(:,:) = dP(:,:) - MMX_I(I,J)*TMP*DM_V(:,:,I)
    enddo
    enddo
    DEALLOCATE(MMX,MMX_I)


    write(6,*) 'maxval(dP)=', maxval(dabs(dP))

    COULOMBV = COULOMBV_SAVE
    DELTAQ = DELTAQ_SAVE
    H0 = H_0
    H = H_SAVE
    FCOUL = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE

    DEALLOCATE(IdentResDM, Res, FCOUL_SAVE, DM_V, UU)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE, BO_SAVE, H_SAVE, H_0)
    DEALLOCATE(X, YY, ORTHOH_SAVE, H_1, H1, dU, nDelta_DO, Delta_DS)
  END SUBROUTINE dP2
!



  SUBROUTINE QMIXPRG(PITER)
#ifdef PROGRESSON
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    INTEGER :: I,J,NUMORB, INDEX
    CHARACTER(20) :: MYMIXERTYPE

    MYMIXERTYPE = MX%MIXERTYPE

    IF(VERBOSE >= 1) WRITE(*,*)"MixerType=", MYMIXERTYPE

    IF(MYMIXERTYPE == "PulayLinear" .AND. PITER >= 10) MYMIXERTYPE = "Linear"

    IF(MYMIXERTYPE == "Linear")THEN

       CALL PRG_LINEARMIXER(DELTAQ,OLDDELTAQS,SCFERROR,MX%MIXCOEFF,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayDM")THEN
       allocate(QTMP1(HDIM*HDIM),QTMP2(HDIM*HDIM))
       QTMP1 = reshape(DOrth,(/HDIM*HDIM/))
       QTMP2 = reshape(DOrth_old,(/HDIM*HDIM/))

       CALL PRG_QMIXER(QTMP1,QTMP2,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
       DOrth = reshape(QTMP1,(/Hdim,Hdim/))
       DOrth_old = reshape(QTMP2,(/Hdim,Hdim/))
       deallocate(QTMP1,QTMP2)

       CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, DOrth, DO_BML, ZERO, LT%MDIM)
       CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, DOrth_old, DO_BML_OLD, ZERO, LT%MDIM)
       CALL BML_COPY_NEW(DO_BML, ORTHOBO_BML)

    ELSEIF(MYMIXERTYPE == "Pulay")THEN

       CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayLinear")THEN

       CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayQlist")THEN

       IF(PITER == 1) OLDQLIST = QLIST
       CALL PRG_QMIXER(QLIST,OLDQLIST,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
       IF(.NOT. ALLOCATED(MYCHARGE)) ALLOCATE(MYCHARGE(NATS))
       INDEX = 0
       MYCHARGE = 0.0d0

       DO I = 1, NATS

          SELECT CASE(BASIS(ELEMPOINTER(I)))

          CASE("s")

             NUMORB = 1

          CASE("p")

             NUMORB = 3

          CASE("d")

             NUMORB = 5

          CASE("f")

             NUMORB = 7

          CASE("sp")

             NUMORB = 4

          CASE("sd")

             NUMORB = 6

          CASE("sf")

             NUMORB = 8

          CASE("pd")

             NUMORB = 8

          CASE("pf")

             NUMORB = 10

          CASE("df")

             NUMORB = 12

          CASE("spd")

             NUMORB = 9

          CASE("spf")

             NUMORB = 11

          CASE("sdf")

             NUMORB = 13

          CASE("pdf")

             NUMORB = 15

          CASE("spdf")

             NUMORB = 16

          END SELECT

          !     MYCHARGE = ZERO
          DO J = 1, NUMORB

             INDEX = INDEX + 1
             MYCHARGE(I) = MYCHARGE(I) + QLIST(INDEX)

          ENDDO

          DELTAQ(I) = MYCHARGE(I) - ATOCC(ELEMPOINTER(I))

       ENDDO

       OLDDELTAQS = DELTAQ

    ELSE
       CALL ERRORS("mixer_mod:qmixprg","Mixing scheme not implemented. &
            & Check MixerType keyword in the input file")
    ENDIF


#endif
  END SUBROUTINE QMIXPRG

END MODULE MIXER_MOD
