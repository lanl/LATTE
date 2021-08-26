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
  USE SPINARRAY
  USE DIAGARRAY    ! CHANGE ANDERS_CHANGE
  USE XBOARRAY     ! CHANGE ANDERS_CHANGE
  USE DMARRAY     ! CHANGE ANDERS
  USE TIMER_MOD
  USE OPENFILES_MOD
#ifdef PROGRESSON
  USE LATTEPARSER
  USE GENXPROGRESS, ONLY : ZMAT_BML, OVER_BML
  USE SPARSEARRAY
  USE PRG_PULAYMIXER_MOD
  USE BML
  USE NONOARRAYPROGRESS
  USE PRG_RESPONSE_MOD
  USE PRG_CHARGES_MOD
#endif
  PRIVATE

#ifdef PROGRESSON
  PUBLIC :: QMIXPRG
#endif

  INTEGER, PARAMETER :: DP = LATTEPREC
  !PUBLIC :: KERNELMIXER  ! this mixer is not used
  PUBLIC :: ADAPTKERNELMIXER
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
    INTEGER :: NDIM, IS
    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:,:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTASPIN_SAVE(:), SUMSPIN_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:,:), H_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:,:)
    !
    NDIM = NATS
    IF (SPINON==1) NDIM = DELTADIM
    ALLOCATE(Res(NDIM*NSPIN),dr(NDIM),vi(NDIM,NDIM),wi(NDIM,NDIM),ui(NDIM,NDIM))
    ALLOCATE(du(NDIM), dq_dv(NDIM,NSPIN), dq_v(NDIM), v(NDIM), ri(NDIM,NDIM))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM,NSPIN), H_SAVE(HDIM,HDIM,2*(NSPIN-1)+1))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(DELTASPIN_SAVE(DELTADIM), SUMSPIN_SAVE(DELTADIM))
    !
#ifdef PROGRESSON
    CALL BML_EXPORT_TO_DENSE(BO_BML,BO)
#elif defined(PROGRESSOFF)
#endif
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    IF (SPINON==1) THEN
      BO_SAVE(:,:,1) = RHOUP
      BO_SAVE(:,:,2) = RHODOWN
      H_SAVE(:,:,1)  = HUP
      H_SAVE(:,:,2)  = HDOWN
      H_SAVE(:,:,3)  = H
      ORTHOH_SAVE(:,:,1) = ORTHOHUP
      ORTHOH_SAVE(:,:,2) = ORTHOHDOWN
      DELTASPIN_SAVE = DELTASPIN
      SUMSPIN_SAVE = SUMSPIN
    ELSE
      BO_SAVE(:,:,1)     = BO
      ORTHOH_SAVE(:,:,1) = ORTHOH
      H_SAVE(:,:,1)      = H
    ENDIF
    H_0 = H0
    H0 = 0.D0

    IF (SPINON==0) THEN
      Res = DELTAQ - PNK(1,:)
      write(6,*) 'MDITER', MDITER, norm2(Res) / SQRT(DBLE(NATS))
    ELSE
      !sumspin = up + down
      !deltaspin = up - down ==>
      ! up = (sum+delta)/2.0
      du(1:NDIM) = DELTASPIN(1:NDIM) - SPIN_PNK(1,1:NDIM)
      dr(1:NDIM) = SUMSPIN(1:NDIM)   - SPIN_PNK(1,1+NDIM:2*NDIM)

      Res(1:NDIM)        = (dr(1:NDIM) + du(1:NDIM)) / 2.D0 ! up
      Res(1+NDIM:2*NDIM) = (dr(1:NDIM) - du(1:NDIM)) / 2.D0 ! down
      write(6,*) 'MDITER', MDITER, norm2(Res(1:NDIM)) / SQRT(DBLE(NDIM)), &
           norm2(Res(1+NDIM:2*NDIM)) / SQRT(DBLE(NDIM))
    ENDIF

    FULL_K = 0.D0

    DO IS = 1, NSPIN
      dr = 0.D0
      do I = 1,NDIM !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
        IF(VERBOSE >= 2)WRITE(*,*)"Constructing response to atom", I,"in FULLKERNELPROPAGATION"
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

        IF (SPINON==0) THEN
          DELTAQ = dq_v
        ELSE
          ! deltapsin is positive/negative for spin up/down
          DELTASPIN = dq_v *(1-2*(IS-1))
          !! get deltaq from deltaspin
          CALL REDUCE_DELTASPIN(NATS,DELTADIM,dq_v,DELTAQ,1)
        ENDIF

        call coulombrspace
        call coulombewald
        CALL ADDQDEP
        IF (SPINON==1) CALL BLDSPINH

        call orthomyh
        Nocc = BNDFIL*float(HDIM)
        beta = 1.D0/KBT
        !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
        IF (NSPIN==1) THEN
          call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
          BO = 2.D0*BO
        ELSE
          call Canon_DM_PRT_SPIN(ORTHOHUP,ORTHOHDOWN,beta,UPEVECS, DOWNEVECS, UPEVALS, DOWNEVALS,CHEMPOT,16,HDIM)
        ENDIF

        call deorthomyrho
        call getdeltaq_resp
        IF (SPINON==1) call getdeltaspin_resp

        IF (SPINON==0) THEN
          dq_dv(:,IS) = DELTAQ
          dr = dq_dv(:,IS)

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
        ELSE
          dq_dv = DELTA_QS
          FULL_K(1     :  NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,1)
          FULL_K(NDIM+1:2*NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,2)
        ENDIF
        !dr(I) = 0.D0
        dr = 0.D0
      enddo
    ENDDO
    IF (SPINON==0) THEN
      FULL_K = 0.D0*FULL_K
      do I = 1,NATS
        !FULL_K(I,I) = MDMIX
        FULL_K(I,I) = 1.0D0
      enddo

      call DGEMM('N','T',NATS,NATS,NATS,1.D0,ui,NATS,wi,NATS,1.D0,FULL_K,NATS)
      ! update q corresponding to q = q - MATMUL(KK,Res)
      !DELTAQ = OLDDELTAQS + QMIX*Res

      dn2dt2(:,1) = MATMUL(FULL_K,Res)
      !!        write(*,*) ' dn2dt2 FULL_K = ', dn2dt2(1:3)
      !!        dn2dt2 = MDMIX*Res
      !!        do I = 1,NATS  !! Let the approximate kernel act on the residual by individual rank-1 updates
      !!           !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
      !!           dn2dt2 = dn2dt2 + dot_product(wi(:,I),Res)*ui(:,I)
      !!        enddo
      !!        write(*,*) ' dn2dt2 Rank-Nats = ', dn2dt2(1:3)
      !!        write(*,*) ' ------------------ '
      !    endif
    ELSE
      FULL_K0 = - FULL_K
      do I = 1,2*NDIM
        FULL_K0(I,I) = FULL_K0(I,I) + 1.0D0
      enddo

      CALL Invert(FULL_K0, FULL_K, 2*NDIM)
      RES = MATMUL(FULL_K,Res)
      dn2dt2(1:NDIM,1) = (Res(1:NDIM) + RES(1+NDIM:2*NDIM))  !sum
      dn2dt2(1:NDIM,2) = (Res(1:NDIM) - RES(1+NDIM:2*NDIM))  !delta
    ENDIF

    IF (SPINON==1) THEN
      RHOUP   = BO_SAVE(:,:,1)
      RHODOWN = BO_SAVE(:,:,2)
      HUP     = H_SAVE(:,:,1)
      HDOWN   = H_SAVE(:,:,2)
      H = H_SAVE(:,:,3)
      ORTHOHUP   = ORTHOH_SAVE(:,:,1)
      ORTHOHDOWN = ORTHOH_SAVE(:,:,2)
      DELTASPIN = DELTASPIN_SAVE
      SUMSPIN = SUMSPIN_SAVE
    ELSE
      BO = BO_SAVE(:,:,1)
      H  = H_SAVE(:,:,1)
      ORTHOH = ORTHOH_SAVE(:,:,1)
    ENDIF
    H0 = H_0
    FCOUL = FCOUL_SAVE
    DELTAQ = DELTAQ_SAVE
    COULOMBV = COULOMBV_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,SUMSPIN_SAVE,DELTASPIN_SAVE)

  END SUBROUTINE FULLKERNELPROPAGATION


#ifdef PROGRESSON
  SUBROUTINE PARALLELFULLKERNELPROPAGATION(MDITER)
    INTEGER, INTENT(IN) :: MDITER
    INTEGER :: I, J, N, ii, jj
    INTEGER :: NDIM, IS, MYIO
    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:,:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTASPIN_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:,:), H_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:,:)
    INTEGER, ALLOCATABLE  :: dummy_array(:), mlsi, mlsii
    REAL(LATTEPREC), ALLOCATABLE :: NUMEL(:)


    TYPE(BML_MATRIX_T) :: zq_bml, zqt_bml, ptham_bml, ptrho_bml, ptaux_bml

    NDIM = NATS

    IF (SPINON==1) NDIM = DELTADIM
    ALLOCATE(Res(NDIM*NSPIN),dr(NDIM),vi(NDIM,NDIM),wi(NDIM,NDIM),ui(NDIM,NDIM))
    ALLOCATE(du(NDIM), dq_dv(NDIM,NSPIN), dq_v(NDIM), v(NDIM), ri(NDIM,NDIM))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM,NSPIN), H_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(DELTASPIN_SAVE(DELTADIM))
    ALLOCATE(NUMEL(NATS))
    ALLOCATE(DUMMY_ARRAY(NATS))
    NUMEL = 0.0d0
    DUMMY_ARRAY = 1
   
    CALL BML_EXPORT_TO_DENSE(BO_BML,BO)

    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    IF (SPINON==1) THEN
      BO_SAVE(:,:,1) = RHOUP
      BO_SAVE(:,:,2) = RHODOWN
      H_SAVE(:,:,1)  = HUP
      H_SAVE(:,:,2)  = HDOWN
      H_SAVE(:,:,3)  = H
      ORTHOH_SAVE(:,:,1) = ORTHOHUP
      ORTHOH_SAVE(:,:,2) = ORTHOHDOWN
      DELTASPIN_SAVE = DELTASPIN
    ELSE
      BO_SAVE(:,:,1)     = BO
      ORTHOH_SAVE(:,:,1) = ORTHOH
      H_SAVE(:,:,1)      = H
    ENDIF
    H_0 = H0
    H0 = 0.D0

    IF (SPINON==0) THEN
      Res = DELTAQ - PNK(1,:)
      write(6,*) 'MDITER', MDITER, norm2(Res) / SQRT(DBLE(NATS))
    ELSE
      !sumspin = up + down
      !deltaspin = up - down ==>
      ! up = (sum+delta)/2.0
      du(1:NDIM) = DELTASPIN(1:NDIM) - SPIN_PNK(1,1:NDIM)
      dr(1:NDIM) = SUMSPIN(1:NDIM)   - SPIN_PNK(1,1+NDIM:2*NDIM)

      Res(1:NDIM)        = (dr(1:NDIM) + du(1:NDIM)) / 2.D0 ! up
      Res(1+NDIM:2*NDIM) = (dr(1:NDIM) - du(1:NDIM)) / 2.D0 ! down
      write(6,*) 'MDITER', MDITER, norm2(Res(1:NDIM)) / SQRT(DBLE(NDIM)), &
           norm2(Res(1+NDIM:2*NDIM)) / SQRT(DBLE(NDIM))
    ENDIF

    FULL_K = 0.D0

    call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptham_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptrho_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,zq_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,zqt_bml)
    call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptaux_bml)
    call bml_multiply(zmat_bml,evecs_bml,zq_bml, 1.0d0,0.0d0,NUMTHRESH)
    call bml_transpose(zq_bml,zqt_bml)

    DO IS = 1, NSPIN
      dr = 0.D0
      do I = 1,NDIM !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
        IF(VERBOSE >= 2)WRITE(*,*)"Constructing response to atom", I,"in PARALLELFULLKERNELPROPAGATION"
        dr(I) = 1.D0
!!!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v
        dq_dv = ZERO
        !dq_v = v/norm2(v)
        dq_v = dr

        IF (SPINON==0) THEN
          DELTAQ = dq_v
        ELSE
          ! deltapsin is positive/negative for spin up/down
          DELTASPIN = dq_v *(1-2*(IS-1))
          !! get deltaq from deltaspin
          CALL REDUCE_DELTASPIN(NATS,DELTADIM,dq_v,DELTAQ,1)
        ENDIF
!        mlsi = time_mls()
        call coulombrspace
        call coulombewald
        CALL ADDQDEP
        call bml_import_from_dense(LT%bml_type, H, ham_bml, NUMTHRESH, HDIM)
        IF (SPINON==1) CALL BLDSPINH
!        write(*,*)"Time for coulombrspace coulombewald ADDQDEP",time_mls() - mlsi
!        mlsi = time_mls()

        call bml_multiply(zqt_bml,ham_bml,ptaux_bml,1.0_dp,0.0_dp,NUMTHRESH)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,NUMTHRESH)
!        write(*,*)"Time for ortho eig",time_mls() - mlsi
!        mlsi = time_mls()
        Nocc = BNDFIL*float(HDIM)
        beta = 1.D0/KBT
        !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
        IF (NSPIN==1) THEN
          !  call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
          !  BO = 2.D0*BO
          call prg_canon_response_orig(ptrho_bml,ptham_bml,nocc,beta,&
               &evals,CHEMPOT,16,numthresh,HDIM)
        ELSE
          call Canon_DM_PRT_SPIN(ORTHOHUP,ORTHOHDOWN,beta,UPEVECS, DOWNEVECS, UPEVALS, DOWNEVALS,CHEMPOT,16,HDIM)
        ENDIF
!        write(*,*)"Time for canon",time_mls() - mlsi

!        mlsi = time_mls()
        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)
!        write(*,*)"Time for deortho eig",time_mls() - mlsi
!        mlsi = time_mls()

        !call bml_export_to_dense(ptrho_bml,BO)
        DELTAQ = 0.0d0
        call prg_get_charges(ptrho_bml, over_bml, NORBINDEX, DELTAQ, numel,&
             &dummy_array, mdim, numthresh)
!        write(*,*)"Time for getcharges",time_mls() - mlsi

        IF (SPINON==1) call getdeltaspin_resp

        IF (SPINON==0) THEN
          dq_dv(:,IS) = DELTAQ
          FULL_K(1:NDIM,I) = dq_dv(1:NDIM,1)
        ELSE
          dq_dv = DELTA_QS
          FULL_K(1     :  NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,1)
          FULL_K(NDIM+1:2*NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,2)
        ENDIF
        dr(I) = 0.D0
      enddo
    ENDDO
    call bml_deallocate(ptrho_bml)
    call bml_deallocate(zq_bml)
    call bml_deallocate(zqt_bml)
    call bml_deallocate(ptaux_bml)
    call bml_deallocate(ptham_bml)
    IF (SPINON==0) THEN
      FULL_K0 = - FULL_K
      do I = 1,NATS
        FULL_K0(I,I) = FULL_K0(I,I) + 1.0D0
      enddo
      CALL Invert(FULL_K0, FULL_K, NATS)

!      mlsii = time_mls()
      dn2dt2(:,1) = MATMUL(FULL_K,Res)
!      write(*,*)"Time for MATMUL",time_mls() - mlsii

    ELSE
      FULL_K0 = - FULL_K
      do I = 1,2*NDIM
        FULL_K0(I,I) = FULL_K0(I,I) + 1.0D0
      enddo

      CALL Invert(FULL_K0, FULL_K, 2*NDIM)
      RES = MATMUL(FULL_K,Res)
      dn2dt2(1:NDIM,1) = (Res(1:NDIM) + RES(1+NDIM:2*NDIM))  !sum
      dn2dt2(1:NDIM,2) = (Res(1:NDIM) - RES(1+NDIM:2*NDIM))  !delta
    ENDIF

    IF (SPINON==1) THEN
      RHOUP   = BO_SAVE(:,:,1)
      RHODOWN = BO_SAVE(:,:,2)
      HUP     = H_SAVE(:,:,1)
      HDOWN   = H_SAVE(:,:,2)
      H = H_SAVE(:,:,3)
      ORTHOHUP   = ORTHOH_SAVE(:,:,1)
      ORTHOHDOWN = ORTHOH_SAVE(:,:,2)
      DELTASPIN = DELTASPIN_SAVE
    ELSE
      BO = BO_SAVE(:,:,1)
      H  = H_SAVE(:,:,1)
      ORTHOH = ORTHOH_SAVE(:,:,1)
    ENDIF
    H0 = H_0
    FCOUL = FCOUL_SAVE
    DELTAQ = DELTAQ_SAVE
    COULOMBV = COULOMBV_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,DELTASPIN_SAVE)
    DEALLOCATE(NUMEL,DUMMY_ARRAY)

  END SUBROUTINE PARALLELFULLKERNELPROPAGATION
#endif


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
    REAL(LATTEPREC), ALLOCATABLE :: FO(:,:), FM(:,:), RI_T(:,:), WORK(:), mlsi
    INTEGER,         ALLOCATABLE :: IPIV(:)
    !
    ALLOCATE(Res(NATS), dr(NATS), vi(NATS,NATS), wi(NATS,NATS), ui(NATS,NATS))
    ALLOCATE(du(NATS), dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(FO(LL,LL), FM(LL,LL), ri_t(LL,NATS),WORK(LL+LL*LL),IPIV(LL))

    !
#ifdef PROGRESSON
    CALL BML_EXPORT_TO_DENSE(BO_BML,BO)
#elif defined(PROGRESSOFF)
#endif
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
      mlsi = time_mls()
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
    dn2dt2(:,1) = 0.D0*Res
    do I = 1,LL
      do J = 1,LL
        dn2dt2(:,1) = dn2dt2(:,1) - FM(I,J)*dot_product(ri(:,J),Res)*vi(:,I)
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
    INTEGER :: NDIM,NDIM2,IS,MYIO

    !REAL(LATTEPREC) :: Res(NATS), dr(NATS), vi(NATS,LL), wi(NATS,LL), ui(NATS,LL)
    !REAL(LATTEPREC) :: dq_dv(NATS), dq_v(NATS), v(NATS), ri(NATS,LL)
    !REAL(LATTEPREC) :: DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS)
    !REAL(LATTEPREC) :: H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM), H_SAVE(HDIM,HDIM)
    REAL(LATTEPREC) :: Nocc, beta, eps
    !REAL(LATTEPREC) :: ORTHOH_SAVE(HDIM,HDIM)
    !REAL(LATTEPREC) :: ri_t(LL,NATS), IDENTRES(NATS)

    REAL(LATTEPREC) :: RESNORM, FEL, DTMP, DIFF, TMP
    INTEGER         :: INFO, RANK
    REAL(LATTEPREC), ALLOCATABLE :: FO(:,:), FM(:,:),WORK(:)
    INTEGER,         ALLOCATABLE :: IPIV(:)

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTASPIN_SAVE(:), SUMSPIN_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:,:), H_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: IDENTRES(:)
    REAL(LATTEPREC) :: mlsi, mls0
    LOGICAL :: COMPUTEKERNEL
    INTEGER, ALLOCATABLE  :: dummy_array(:)
    REAL(LATTEPREC), ALLOCATABLE :: NUMEL(:)

#ifdef PROGRESSON
    TYPE(BML_MATRIX_T) :: zq_bml, zqt_bml, ptham_bml, ptrho_bml, ptaux_bml
#endif

    !
    NDIM = NATS
    IF (SPINON==1) NDIM = DELTADIM
    NDIM2 = NDIM * NSPIN
    !
    ALLOCATE(Res(NDIM2),dr(NDIM2),vi(NDIM2,NDIM2))
    ALLOCATE(dq_dv(NDIM2), dq_v(NDIM), v(NDIM2), ri(NDIM2,NDIM2))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM,NSPIN), H_SAVE(HDIM,HDIM,2*(NSPIN-1)+1))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(WORK(LL+LL*LL),IDENTRES(NDIM2))
    ALLOCATE(DELTASPIN_SAVE(DELTADIM),SUMSPIN_SAVE(DELTADIM))
    !
#ifdef PROGRESSON
    CALL BML_EXPORT_TO_DENSE(BO_BML,BO)
#elif defined(PROGRESSOFF)
#endif

    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL

    IF (SPINON==1) THEN
      !SAVE RHOUP/DOWN, H_UP/DOWN
      !ORTHOHUP/DOWN
      ORTHOH_SAVE(:,:,1) = ORTHOHUP
      ORTHOH_SAVE(:,:,2) = ORTHOHDOWN
      H_SAVE(:,:,1) = HUP
      H_SAVE(:,:,2) = HDOWN
      H_SAVE(:,:,3) = H
      BO_SAVE(:,:,1) = RHOUP
      BO_SAVE(:,:,2) = RHODOWN
      DELTASPIN_SAVE = DELTASPIN
      SUMSPIN_SAVE = SUMSPIN
    ELSE
      BO_SAVE(:,:,1)     = BO
      ORTHOH_SAVE(:,:,1) = ORTHOH
      H_SAVE(:,:,1)      = H
    ENDIF
    H_0 = H0
    H0 = 0.D0

    IF (SPINON==0) THEN
      Res = DELTAQ - PNK(1,:)
      RESNORM = NORM2(RES) / SQRT(DBLE(NATS))
      write(6,*) 'MDITER, RESNORM', MDITER, RESNORM
    ELSE
      dr(1:NDIM)        = SUMSPIN(1:NDIM)   - SPIN_PNK(1,1+NDIM:2*NDIM)
      dr(1+NDIM:2*NDIM) = DELTASPIN(1:NDIM) - SPIN_PNK(1,1:NDIM)

      Res(1:NDIM)        = (dr(1:NDIM) + dr(1+NDIM:2*NDIM)) / 2.D0 ! up
      Res(1+NDIM:2*NDIM) = (dr(1:NDIM) - dr(1+NDIM:2*NDIM)) / 2.D0 ! down
      write(6,*) 'MDITER', MDITER, norm2(Res(1:NDIM)) / SQRT(DBLE(NDIM)), &
           norm2(Res(1+NDIM:2*NDIM)) / SQRT(DBLE(NDIM))
    ENDIF

    RANK = 0
    I = 0

    IF (MDITER == 1) THEN

      IF(READKERNEL)THEN
        WRITE(*,*)"Reading the kernel from file ..."
        CALL OPEN_FILE_TO_READ(MYIO,"kernel.tmp")
        DIFF = 0.0d0
        READ(MYIO,*)TMP
        IF(TMP .NE. HDIM)THEN
          WRITE(*,*)"WARNING: The kernel.tmp file is inconsistent with this system &
            & I will recompute the Kernel instead ..."
          COMPUTEKERNEL = .TRUE.
        ELSE 
          DO I=1,HDIM
            READ(MYIO,*)TMP
            DIFF = DIFF + ABS(TMP - H(I,1))
          ENDDO
          IF(DIFF > 1.0E-10)then
            WRITE(*,*)"WARNING: The kernel.tmp file is inconsistent with this system &
              & I will recompute the Kernel instead ..."
            COMPUTEKERNEL = .TRUE.
          ELSE  
            DO I=1,NATS
              DO J=1,NATS
                READ(MYIO,*)FULL_K(I,J)
              ENDDO
            ENDDO
            COMPUTEKERNEL = .FALSE.
          ENDIF
        ENDIF
        CLOSE(MYIO)
      ELSE
        COMPUTEKERNEL = .TRUE.
      ENDIF 

      IF(COMPUTEKERNEL)THEN
#ifdef PROGRESSON
        CALL PARALLELFULLKERNELPROPAGATION(MDITER)
#else
        CALL FULLKERNELPROPAGATION(MDITER)
#endif
        IF(SAVEKERNEL)THEN
          WRITE(*,*)"Saving kernel into file ..."
          CALL OPEN_FILE(MYIO,"kernel.tmp",.TRUE.)
          !Creating a "tag" with the Hamiltonian - for checking purposes
          WRITE(MYIO,*)HDIM
          DO I=1,HDIM
            WRITE(MYIO,*)H(I,1)
          ENDDO
          DO I=1,NATS
            DO J=1,NATS
              WRITE(MYIO,*)FULL_K(I,J)
            ENDDO
          ENDDO
          CLOSE(MYIO)
        ENDIF
      ENDIF

    ELSE
      MLS0 = time_mls()
      Res = MATMUL(FULL_K,Res) !! FULL_KK is the preconditioner
      dr = Res

      I = 0
      FEL = 1.D0

#ifdef PROGRESSON
      call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptham_bml)
      call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptrho_bml)
      call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,zq_bml)
      call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,zqt_bml)
      call bml_zero_matrix(lt%bml_type,bml_element_real,LATTEPREC,HDIM,HDIM,ptaux_bml)
      call bml_multiply(zmat_bml,evecs_bml,zq_bml, 1.0d0,0.0d0,NUMTHRESH)
      call bml_transpose(zq_bml,zqt_bml)
      ALLOCATE(NUMEL(NATS))
      ALLOCATE(DUMMY_ARRAY(NATS))
      NUMEL = 0.0d0
      DUMMY_ARRAY = 1
#endif

      RANK = 0
      DO WHILE ((FEL > KERNELTOL) .AND. (RANK <= LL))  !! LL is the number of rank-1 updates  LL = 0 means preconditioning only!
        WRITE(*,*)"Adapting the Kernel, FEL > KERNELTOL ...",I
        I = I + 1
!        mlsi = time_mls()
        vi(:,I) = dr/norm2(dr)
        do J = 1,I-1
          vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
        enddo
        vi(:,I) = vi(:,I)/norm2(vi(:,I))
        v(:) = vi(:,I)
!!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v

        v = v / norm2(v)
!        write(*,*)"Time for getting vs at rankN",time_mls() - mlsi

        dq_dv = ZERO
        IF (SPINON==1) THEN
          DELTASPIN    = v(1:NDIM)  - v(1+NDIM:2*NDIM)
          dq_v(1:NDIM) = v(1:NDIM)  + v(1+NDIM:2*NDIM)
!!!! get deltaq from deltaspin
          CALL REDUCE_DELTASPIN(NATS,DELTADIM,dq_v,DELTAQ,1)
        ELSE
          dq_v   = v
          DELTAQ = dq_v
        ENDIF

!        mlsi = time_mls()
        call coulombrspace
        call coulombewald
        call addqdep
!        write(*,*)"Time for coulombrspace coulombewald addqdep at rankN",time_mls() - mlsi

        IF (SPINON==1) CALL BLDSPINH

!        mlsi = time_mls()
#ifdef PROGRESSON
        call bml_import_from_dense(LT%bml_type, H, ham_bml, NUMTHRESH, HDIM)
        call bml_multiply(zqt_bml,ham_bml,ptaux_bml,1.0_dp,0.0_dp,NUMTHRESH)
        call bml_multiply(ptaux_bml,zq_bml,ptham_bml,1.0_dp,0.0_dp,NUMTHRESH)
#else
        call orthomyh
#endif        
!        write(*,*)"Time for orthomyh at rankN",time_mls() - mlsi

        Nocc = BNDFIL*float(HDIM)
        beta = 1.D0/KBT

        IF (SPINON==0) THEN

!          mlsi = time_mls()
#ifdef PROGRESSON
          call prg_canon_response_orig(ptrho_bml,ptham_bml,nocc,beta,&
               &evals,CHEMPOT,16,numthresh,HDIM)
#else
          call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
 !         call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
  !        BO = 2.D0*BO
#endif             
!          write(*,*)"Time for Canon_DM_PRT at rankN",time_mls() - mlsi
        ELSE
          call Canon_DM_PRT_SPIN(ORTHOHUP,ORTHOHDOWN,beta,UPEVECS, DOWNEVECS, UPEVALS, DOWNEVALS,CHEMPOT,16,HDIM)
        ENDIF

!        mlsi = time_mls()
#ifdef PROGRESSON
        call bml_multiply(zq_bml,ptrho_bml,ptaux_bml,1.0_dp,0.0_dp,0.0_dp)
        call bml_multiply(ptaux_bml,zqt_bml,ptrho_bml,2.0_dp,0.0_dp,0.0_dp)
#else
        call deorthomyrho
#endif             
!        write(*,*)"Time for deorthomyrho at rankN",time_mls() - mlsi

!        mlsi = time_mls()
#ifdef PROGRESSON
        call prg_get_charges(ptrho_bml, over_bml, NORBINDEX, DELTAQ, numel,&
             &dummy_array, hdim, numthresh)
#else
        call getdeltaq_resp
#endif             
!        write(*,*)"Time for getdeltaq_resp at rankN",time_mls() - mlsi
        
        IF (SPINON==1) call getdeltaspin_resp

        IF (SPINON==0) THEN
          dq_dv = DELTAQ
        ELSE
          dq_dv(1     :  NDIM) = DELTA_QS(1:NDIM,1)
          dq_dv(1+NDIM:2*NDIM) = DELTA_QS(1:NDIM,2)
        ENDIF

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
        v        = 0.D0
        DO K = 1,RANK
          DO J = 1,RANK
            DTMP = FM(K,J)*dot_product(RI(:,J),RES)
            IdentRes = IdentRes + DTMP*RI(:,K)
            v(:) = v(:) - DTMP*VI(:,K)
          ENDDO
        ENDDO

        IF (SPINON==0) THEN
          dn2dt2(:,1) = v(:)
        ELSE
          dn2dt2(1:NDIM,1) = (v(1:NDIM) + v(1+NDIM:2*NDIM)) !/ 2.D0
          dn2dt2(1:NDIM,2) = (v(1:NDIM) - v(1+NDIM:2*NDIM)) !/ 2.D0
        ENDIF

        FEL = NORM2(IDENTRES - RES) / NORM2(RES)
        !WRITE(6,*) 'ERROR:', MDITER, I, FEL

        DEALLOCATE(FO, FM, IPIV)

      ENDDO

#ifdef PROGRESSON
      call bml_deallocate(ptrho_bml)
      call bml_deallocate(zq_bml)
      call bml_deallocate(zqt_bml)
      call bml_deallocate(ptaux_bml)
      call bml_deallocate(ptham_bml)
#endif

    COULOMBV = COULOMBV_SAVE
    IF (SPINON==1) THEN
      RHOUP   = BO_SAVE(:,:,1)
      RHODOWN = BO_SAVE(:,:,2)
      HUP     = H_SAVE(:,:,1)
      HDOWN   = H_SAVE(:,:,2)
      H = H_SAVE(:,:,3)
      ORTHOHUP   = ORTHOH_SAVE(:,:,1)
      ORTHOHDOWN = ORTHOH_SAVE(:,:,2)
      DELTASPIN = DELTASPIN_SAVE
      SUMSPIN = SUMSPIN_SAVE
    ELSE
      H      = H_SAVE(:,:,1)
      BO     = BO_SAVE(:,:,1)
      ORTHOH = ORTHOH_SAVE(:,:,1)
    ENDIF
    H0 = H_0
    FCOUL = FCOUL_SAVE
    DELTAQ = DELTAQ_SAVE


    DEALLOCATE(RES, DR, VI, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,WORK,IDENTRES,DELTASPIN_SAVE, SUMSPIN_SAVE)

#ifdef PROGRESSON      
    IF (ALLOCATED(NUMEL)) DEALLOCATE(NUMEL,DUMMY_ARRAY)
#endif      

    IF (RANK >= LL ) THEN
#ifdef PROGRESSON
       CALL PARALLELFULLKERNELPROPAGATION(MDITER)
#else
       CALL FULLKERNELPROPAGATION(MDITER)
#endif
    ENDIF


      write(*,*)"Time for rankN update", time_mls() - mls0,"With iter",I
    endif


  END SUBROUTINE ADAPTPRECONDKERNEL

  SUBROUTINE FULLKERNELMIXER(NDIM, dspin)
    INTEGER, INTENT(IN) :: NDIM
    REAL(LATTEPREC), INTENT(OUT) :: dspin(NDIM, NSPIN)

    INTEGER :: I, J, N, ii, jj
    INTEGER :: IS
    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:), wi(:,:), ui(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: du(:), dq_dv(:,:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTASPIN_SAVE(:), SUMSPIN_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:,:), H_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:,:)
    !
    ALLOCATE(Res(NDIM*NSPIN),dr(NDIM),vi(NDIM,NDIM),wi(NDIM,NDIM),ui(NDIM,NDIM))
    ALLOCATE(du(NDIM), dq_dv(NDIM,NSPIN), dq_v(NDIM), v(NDIM), ri(NDIM,NDIM))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM,NSPIN), H_SAVE(HDIM,HDIM,2*(NSPIN-1)+1))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(DELTASPIN_SAVE(DELTADIM), SUMSPIN_SAVE(DELTADIM))
    !
#ifdef PROGRESSON
     CALL BML_EXPORT_TO_DENSE(BO_BML,BO)
#elif defined(PROGRESSOFF)
#endif
    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL
    IF (SPINON==1) THEN
      BO_SAVE(:,:,1) = RHOUP
      BO_SAVE(:,:,2) = RHODOWN
      H_SAVE(:,:,1)  = HUP
      H_SAVE(:,:,2)  = HDOWN
      H_SAVE(:,:,3)  = H
      ORTHOH_SAVE(:,:,1) = ORTHOHUP
      ORTHOH_SAVE(:,:,2) = ORTHOHDOWN
      DELTASPIN_SAVE = DELTASPIN
      SUMSPIN_SAVE = SUMSPIN
    ELSE
      BO_SAVE(:,:,1)     = BO
      ORTHOH_SAVE(:,:,1) = ORTHOH
      H_SAVE(:,:,1)      = H
    ENDIF
    H_0 = H0
    H0 = 0.D0

    IF (SPINON==0) THEN
      Res = DELTAQ - OLDDELTAQS
      write(6,*) 'Residue=', norm2(Res) / SQRT(DBLE(NATS)) 
    ELSE
      !sumspin = up + down
      !deltaspin = up - down ==>
      ! up = (sum+delta)/2.0
      du(1:NDIM) = DELTASPIN(1:NDIM) - OLDDELTASPIN(1:NDIM)
      dr(1:NDIM) = SUMSPIN(1:NDIM) - OLDSUMSPIN(1:NDIM)

      Res(1:NDIM)        = (dr(1:NDIM) + du(1:NDIM)) / 2.D0 ! up
      Res(1+NDIM:2*NDIM) = (dr(1:NDIM) - du(1:NDIM)) / 2.D0 ! down
      write(6,*) 'Residue=', norm2(Res(1:2*NDIM)) / SQRT(DBLE(2*NDIM))

      write(6,*) 'Residue=', norm2(Res(1:NDIM)) / SQRT(DBLE(NDIM)), &
                            norm2(Res(1+NDIM:2*NDIM)) / SQRT(DBLE(NDIM)) 
    ENDIF

    FULL_S = 0.D0

    DO I = 1,NDIM !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
    !DO IS = 1, NSPIN
        dr = 0.D0
!! OMP not done yet!
!!$OMP PARALLEL DO DEFAULT (NONE) &
!!$OMP PRIVATE(I,dq_dv, DELTAQ, NOCC,BETA) 
        !do I = 1,NDIM !! NATS is the number of rank-1 updates  LL = 0 means linear mixing
        DO IS = 1, NSPIN
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

           IF (SPINON==0) THEN
             DELTAQ = dq_v
           ELSE
             ! deltapsin is positive/negative for spin up/down
             DELTASPIN = dq_v *(1-2*(IS-1))
             !! get deltaq from deltaspin
             CALL REDUCE_DELTASPIN(NATS,DELTADIM,dq_v,DELTAQ,1)
           ENDIF

           !IF (IS==1) call coulombrspace
           !IF (IS==1) call coulombewald
           call coulombrspace
           call coulombewald

           CALL ADDQDEP
           IF (SPINON==1) CALL BLDSPINH

           call orthomyh
           Nocc = BNDFIL*float(HDIM)
           beta = 1.D0/KBT
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           IF (NSPIN==1) THEN
             call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)
             BO = 2.D0*BO
           ELSE
             call Canon_DM_PRT_SPIN(ORTHOHUP,ORTHOHDOWN,beta,UPEVECS, DOWNEVECS, UPEVALS, DOWNEVALS,CHEMPOT,16,HDIM)
           ENDIF

           call deorthomyrho
           call getdeltaq_resp
           IF (SPINON==1) call getdeltaspin_resp

           IF (SPINON==0) THEN
             dq_dv(:,IS) = DELTAQ
             dr = dq_dv(:,IS)

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
           ELSE
             dq_dv = DELTA_QS 
             FULL_S(1     :  NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,1)
             FULL_S(NDIM+1:2*NDIM,I+NDIM*(IS-1)) = dq_dv(1:NDIM,2)

             !FULL_S(1+NDIM*(IS-1):NDIM*IS, I+NDIM*(IS-1)) = dq_dv(1:NDIM,IS)
           ENDIF
           !dr(I) = 0.D0
           dr = 0.D0
        enddo
    ENDDO
!!$OMP END PARALLEL DO

    IF (SPINON==0) THEN
        FULL_S = 0.D0
        do I = 1,NATS 
           FULL_S(I,I) = 1.0D0
        enddo

        call DGEMM('N','T',NATS,NATS,NATS,1.D0,ui,NATS,wi,NATS,1.D0,FULL_S,NATS)
        ! update q corresponding to q = q - MATMUL(KK,Res)
        !DELTAQ = OLDDELTAQS + QMIX*Res

        dspin(:,1) = MATMUL(FULL_S,Res)
    ELSE
      FULL_K0 = - FULL_S
      do I = 1,2*NDIM 
         FULL_K0(I,I) = FULL_K0(I,I) + 1.0D0
      enddo
      
      CALL Invert(FULL_K0, FULL_S, 2*NDIM)
      RES = MATMUL(FULL_S, Res)
      dspin(1:NDIM,1) = (Res(1:NDIM) + RES(1+NDIM:2*NDIM)) !/ 2.D0 !sum
      dspin(1:NDIM,2) = (Res(1:NDIM) - RES(1+NDIM:2*NDIM)) !/ 2.D0 !delta
    ENDIF

    IF (SPINON==1) THEN
      RHOUP   = BO_SAVE(:,:,1)
      RHODOWN = BO_SAVE(:,:,2)
      HUP     = H_SAVE(:,:,1)
      HDOWN   = H_SAVE(:,:,2)
      H = H_SAVE(:,:,3)
      ORTHOHUP   = ORTHOH_SAVE(:,:,1)
      ORTHOHDOWN = ORTHOH_SAVE(:,:,2)
      DELTASPIN = DELTASPIN_SAVE
      SUMSPIN = SUMSPIN_SAVE
    ELSE
      BO = BO_SAVE(:,:,1)
      H  = H_SAVE(:,:,1)
      ORTHOH = ORTHOH_SAVE(:,:,1)
    ENDIF
    H0 = H_0
    FCOUL = FCOUL_SAVE
    DELTAQ = DELTAQ_SAVE
    COULOMBV = COULOMBV_SAVE

    DEALLOCATE(RES, DR, VI, WI, UI, DU, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,SUMSPIN_SAVE,DELTASPIN_SAVE)

  END SUBROUTINE FULLKERNELMIXER

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
      dn2dt2(:,1) = MDMIX*Res
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
      dn2dt2(:,1) = MDMIX*Res
      do I = 1,LL  !! Let the approximate kernel act on the residual by individual rank-1 updates
        !DELTAQ = DELTAQ + dot_product(wi(:,I),Res)*ui(:,I)
        dn2dt2(:,1) = dn2dt2(:,1) + dot_product(wi(:,I),Res)*ui(:,I)
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

  SUBROUTINE ADAPTKERNELMIXER(FIRSTCALL,LL)
    implicit none
    LOGICAL, INTENT(IN) :: FIRSTCALL
    INTEGER, INTENT(IN) :: LL
    INTEGER :: I, J, N, II, JJ, K
    INTEGER :: NDIM,NDIM2,IS

    REAL(LATTEPREC) :: Nocc, beta, eps

    REAL(LATTEPREC) :: RESNORM, FEL, DTMP
    INTEGER         :: INFO, RANK
    REAL(LATTEPREC), ALLOCATABLE :: FO(:,:), FM(:,:),WORK(:)
    INTEGER,         ALLOCATABLE :: IPIV(:)

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), dr(:), vi(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: dq_dv(:), dq_v(:), v(:), ri(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTASPIN_SAVE(:), SUMSPIN_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), BO_SAVE(:,:,:), H_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: FCOUL_SAVE(:,:), ORTHOH_SAVE(:,:,:)
    REAL(LATTEPREC), ALLOCATABLE :: dspin(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: IDENTRES(:)
    !
    NDIM = NATS
    IF (SPINON==1) NDIM = DELTADIM
    NDIM2 = NDIM * NSPIN
    !
    ALLOCATE(Res(NDIM2),dr(NDIM2),vi(NDIM2,NDIM2))
    ALLOCATE(dq_dv(NDIM2), dq_v(NDIM), v(NDIM2), ri(NDIM2,NDIM2))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), BO_SAVE(HDIM,HDIM,NSPIN), H_SAVE(HDIM,HDIM,2*(NSPIN-1)+1))
    ALLOCATE(FCOUL_SAVE(3,NATS),ORTHOH_SAVE(HDIM,HDIM,NSPIN))
    ALLOCATE(WORK(LL+LL*LL),IDENTRES(NDIM2))
    ALLOCATE(DELTASPIN_SAVE(DELTADIM),SUMSPIN_SAVE(DELTADIM))
    ALLOCATE(dspin(NDIM,NSPIN))

    !
#ifdef PROGRESSON
     CALL BML_EXPORT_TO_DENSE(BO_BML,BO)
#elif defined(PROGRESSOFF)
#endif

    DELTAQ_SAVE = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE = FCOUL

    IF (SPINON==1) THEN
      !SAVE RHOUP/DOWN, H_UP/DOWN
      !ORTHOHUP/DOWN
      ORTHOH_SAVE(:,:,1) = ORTHOHUP
      ORTHOH_SAVE(:,:,2) = ORTHOHDOWN
      H_SAVE(:,:,1) = HUP
      H_SAVE(:,:,2) = HDOWN
      H_SAVE(:,:,3) = H
      BO_SAVE(:,:,1) = RHOUP
      BO_SAVE(:,:,2) = RHODOWN
      DELTASPIN_SAVE = DELTASPIN
      SUMSPIN_SAVE = SUMSPIN
    ELSE
      BO_SAVE(:,:,1)     = BO
      ORTHOH_SAVE(:,:,1) = ORTHOH
      H_SAVE(:,:,1)      = H
    ENDIF
    H_0 = H0
    H0 = 0.D0

    !!!
    !!! Feb 26, 2021
    !!! important problem, deltaspin's dim is different deltaq's
    !!! but deltaq can be calculated q_up and q_down
    !!! so, we have to use q_up/down as variables for spin polarized
    !!! calculations
    !!!
    IF (SPINON==0) THEN
      Res = DELTAQ - OLDDELTAQS
      RESNORM = NORM2(RES) / SQRT(DBLE(NATS)) 
      write(6,*) 'ITER', RESNORM
    ELSE
      dr(1:NDIM)        = SUMSPIN(1:NDIM)   - OLDSUMSPIN(1:NDIM)
      dr(1+NDIM:2*NDIM) = DELTASPIN(1:NDIM) - OLDDELTASPIN(1:NDIM)
      
      Res(1:NDIM)        = (dr(1:NDIM) + dr(1+NDIM:2*NDIM)) / 2.D0 ! up
      Res(1+NDIM:2*NDIM) = (dr(1:NDIM) - dr(1+NDIM:2*NDIM)) / 2.D0 ! down
      write(6,*) 'ITER', norm2(Res(1:NDIM)) / SQRT(DBLE(NDIM)), &
                         norm2(Res(1+NDIM:2*NDIM)) / SQRT(DBLE(NDIM)), &
                         SUM(DELTASPIN(1:NDIM)) 
    ENDIF


    IF (FIRSTCALL) THEN
       write(6,*) 'Test: pre-condition for SCF mixer!'
       CALL FULLKERNELMIXER(NDIM, dspin)
    ELSE
       write(6,*) 'Test: adaptkernel mixer!'
      Res = MATMUL(FULL_S, Res) !! FULL_S is the preconditioner
      dr = Res

      I = 0
      RANK = 0
      FEL = 1.D0
      DO WHILE (FEL > KERNELTOL .AND. RANK <= LL)  !! LL is the number of rank-1 updates  LL = 0 means preconditioning only!
         I = I + 1

         ! use first half of dr for spin up and the second half for down???
         vi(:,I) = dr/norm2(dr)
         do J = 1,I-1
            vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)
         enddo
         vi(:,I) = vi(:,I)/norm2(vi(:,I))
         v(:) = vi(:,I)
         !!!! Calculated dq_dv, which is the response in q(n) from change in input charge n = v

         v = v / norm2(v)

         dq_dv = ZERO
         IF (SPINON==1) THEN
           DELTASPIN    = v(1:NDIM)  - v(1+NDIM:2*NDIM)
           dq_v(1:NDIM) = v(1:NDIM)  + v(1+NDIM:2*NDIM)
           !!!! get deltaq from deltaspin
           CALL REDUCE_DELTASPIN(NATS,DELTADIM,dq_v,DELTAQ,1)
         ELSE
           dq_v   = v
           DELTAQ = dq_v
         ENDIF

         call coulombrspace
         call coulombewald
         call addqdep
         IF (SPINON==1) CALL BLDSPINH
         call orthomyh
         Nocc = BNDFIL*float(HDIM)
         beta = 1.D0/KBT

         IF (SPINON==0) THEN
           !call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
           call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,18,HDIM)
           BO = 2.D0*BO
         ELSE
           call Canon_DM_PRT_SPIN(ORTHOHUP,ORTHOHDOWN,beta,UPEVECS, DOWNEVECS, UPEVALS, DOWNEVALS,CHEMPOT,18,HDIM)
         ENDIF

         call deorthomyrho
         call getdeltaq_resp
         IF (SPINON==1) call getdeltaspin_resp

         IF (SPINON==0) THEN
           dq_dv = DELTAQ
         ELSE
           dq_dv(1     :  NDIM) = DELTA_QS(1:NDIM,1)
           dq_dv(1+NDIM:2*NDIM) = DELTA_QS(1:NDIM,2)
         ENDIF

         dr = dq_dv - v
         dr = MATMUL(FULL_S, dr)
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

         FM       = FO
         DN2DT2   = 0.D0
         IDENTRES = 0.D0
         v        = 0.D0
         DO K = 1,RANK
         DO J = 1,RANK
           DTMP = FM(K,J)*dot_product(RI(:,J),RES)
           IdentRes = IdentRes + DTMP*RI(:,K)
           v(:) = v(:) - DTMP*VI(:,K)
         ENDDO
         ENDDO

         IF (SPINON==0) THEN
            dspin(:,1) = v(:)
         ELSE
            dspin(1:NDIM,1) = (v(1:NDIM) + v(1+NDIM:2*NDIM)) !/ 2.D0
            dspin(1:NDIM,2) = (v(1:NDIM) - v(1+NDIM:2*NDIM)) !/ 2.D0
         ENDIF

         FEL = NORM2(IDENTRES - RES) / NORM2(RES)
         WRITE(6,*) 'ERROR:', I, FEL

         DEALLOCATE(FO, FM, IPIV)

      ENDDO 

    endif

    COULOMBV = COULOMBV_SAVE
    IF (SPINON==1) THEN
      RHOUP   = BO_SAVE(:,:,1)
      RHODOWN = BO_SAVE(:,:,2)
      HUP     = H_SAVE(:,:,1)
      HDOWN   = H_SAVE(:,:,2)
      H = H_SAVE(:,:,3)
      ORTHOHUP   = ORTHOH_SAVE(:,:,1)
      ORTHOHDOWN = ORTHOH_SAVE(:,:,2)
      DELTASPIN = DELTASPIN_SAVE
      SUMSPIN = SUMSPIN_SAVE
    ELSE
      H      = H_SAVE(:,:,1)
      BO     = BO_SAVE(:,:,1)
      ORTHOH = ORTHOH_SAVE(:,:,1)
    ENDIF
    H0 = H_0
    FCOUL = FCOUL_SAVE
    DELTAQ = DELTAQ_SAVE

    IF (RANK >= LL ) THEN
       WRITE(6,*) 'Update kernel for mixer!'
       CALL FULLKERNELMIXER(NDIM, dspin)
    ENDIF
    
    IF (SPINON==1) THEN
      !DELTASPIN = OLDDELTASPIN + QMIX * (Res(1:NDIM) - Res(1+NDIM:2*NDIM))/2.d0
      !SUMSPIN = OLDSUMSPIN + QMIX * (Res(1:NDIM) + Res(1+NDIM:2*NDIM))/2.d0
      !DELTASPIN = DELTASPIN + dspin(1:NDIM,2)
      !SUMSPIN = SUMSPIN + dspin(1:NDIM,1)

      DELTASPIN = OLDDELTASPIN + dspin(1:NDIM,2)
      SUMSPIN = OLDSUMSPIN + dspin(1:NDIM,1)

    ELSE
      DELTAQ = OLDDELTAQS !+ QMIX * Res(1:NDIM)
      DELTAQ = DELTAQ + dspin(1:NDIM,1)
    ENDIF

    DEALLOCATE(RES, DR, VI, DQ_DV, DQ_V, V, RI)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE,H_0, BO_SAVE, H_SAVE, dspin)
    DEALLOCATE(FCOUL_SAVE,ORTHOH_SAVE,WORK,IDENTRES,SUMSPIN_SAVE, DELTASPIN_SAVE)
  END SUBROUTINE ADAPTKERNELMIXER

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
    !
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
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:), Delta_DS(:,:)
#ifdef PROGRESSON
    TYPE(BML_MATRIX_T) :: AUX_BML, X_BML, Y_BML, BO_BML_SAVE
    TYPE(BML_MATRIX_T) :: DELTA_BML
#endif

    ALLOCATE(Res(NATS),BO_SAVE(HDIM,HDIM),X(HDIM,HDIM), YY(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM), Delta_DS(HDIM,HDIM))

    !
    BO_SAVE = BO
#ifdef PROGRESSON
    CALL BML_EXPORT_TO_DENSE(PO0_BML, PO_0)
    CALL BML_EXPORT_TO_DENSE(DO_BML, DOrth)

    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,X_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,Y_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,BO_BML_SAVE)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,DELTA_BML)

    !CALL BML_COPY_NEW(ORTHOBO_BML, BO_BML_SAVE)
    CALL BML_COPY_NEW(DO_BML, DELTA_BML)

    CALL BML_ADD(DELTA_BML,PO0_BML,1.0_DP,-1.0_DP,NUMTHRESH)
    CALL BML_EXPORT_TO_DENSE(DELTA_BML, Delta_DO)

    CALL BML_MULTIPLY(DELTA_BML,ZMAT_BML,Y_BML,1.0_DP,0.0_DP,NUMTHRESH)
    CALL BML_TRANSPOSE(ZMAT_BML, AUX_BML)
    CALL BML_MULTIPLY(Y_BML,AUX_BML, X_BML,  1.0_DP,0.0_DP,NUMTHRESH)
    CALL BML_MULTIPLY(X_BML,OVER_BML,AUX_BML,1.0_DP,0.0_DP,NUMTHRESH)

    CALL BML_EXPORT_TO_DENSE(AUX_BML, Delta_DS)
#elif defined(PROGRESSOFF)
    Delta_DO = DOrth-PO_0

    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,Delta_DO,HDIM,ZERO,YY,HDIM)
    call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,YY,HDIM,XMAT,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,X,HDIM,SMAT,HDIM,ZERO,Delta_DS,HDIM)
#endif

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
#ifdef PROGRESSON
      CALL BML_COPY_NEW(DELTA_BML, D2P_BML)
#elif defined(PROGRESSOFF)
      d2PO = Delta_DO
#endif
      !write(6,*) PO_0
    ELSE
      write(*,*) ' RESPONSE MIX MIXER_MOD (MD)'
#ifdef PROGRESSON
      CALL dP2PRG(DELTA_BML, d2P_BML)
#elif defined(PROGRESSOFF)
      call dP2(Delta_DO, d2PO)
#endif
    ENDIF

    BO = BO_SAVE

    DEALLOCATE(RES, BO_SAVE, X, YY, DELTA_DO, DELTA_DS)

#ifdef PROGRESSON
    !CALL BML_COPY_NEW(BO_BML_SAVE, ORTHOBO_BML)
    CALL BML_DEALLOCATE(BO_BML_SAVE)
    CALL BML_DEALLOCATE(DELTA_BML)
    CALL BML_DEALLOCATE(AUX_BML)
    CALL BML_DEALLOCATE(X_BML)
    CALL BML_DEALLOCATE(Y_BML)
#endif

  END SUBROUTINE dP2MD
  !
  SUBROUTINE dP2MIXER(PITER,SCF_ERR,MAXDQ)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER
    REAL(LATTEPREC), INTENT(IN) :: SCF_ERR
    REAL(LATTEPREC), INTENT(IN) :: MAXDQ
    !
    INTEGER :: I, J, K, N, ii, jj
#ifdef PROGRESSON
    TYPE(BML_MATRIX_T) :: AUX_BML, DU_BML
#endif

    REAL(LATTEPREC), ALLOCATABLE :: DU(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DO(:,:)

    ALLOCATE(DU(HDIM,HDIM))
    ALLOCATE(Delta_DO(HDIM,HDIM))
    !
#ifdef PROGRESSON
    !CALL BML_EXPORT_TO_DENSE(DO_BML, DOrth)
    !CALL BML_EXPORT_TO_DENSE(DO_BML_OLD, DOrth_old)

    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,DU_BML)

    CALL BML_COPY_NEW(DO_BML, AUX_BML)
    CALL BML_COPY_NEW(DO_BML_OLD, ORTHOBO_BML)
    CALL BML_ADD(AUX_BML,DO_BML_OLD,1.0_DP,-1.0_DP,NUMTHRESH)

    CALL BML_EXPORT_TO_DENSE(AUX_BML, Delta_DO)
#elif defined(PROGRESSOFF)
    Delta_DO = DOrth-DOrth_old
#endif

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
#ifdef PROGRESSON
    IF(SCF_ERR > 0.1D0) THEN
      CALL BML_COPY_NEW(AUX_BML, DU_BML)
      CALL BML_ADD(DU_BML, AUX_BML, QMIX, 0.D0, 1.D-20)
      !dU = QMIX*Delta_DO
    ELSE
      !call dP2(Delta_DO, dU)
      call dP2PRG(AUX_BML, DU_BML)
    ENDIF

    !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, DU, DU_BML, ZERO, LT%MDIM)
    CALL BML_ADD(ORTHOBO_BML,DU_BML,1.0_DP,1.0_DP,NUMTHRESH)
#elif defined(PROGRESSOFF)
    !write(*,*) ' LINEAR MIX MIXER_MOD', QMIX, SCF_ERR
    IF(SCF_ERR > 0.1D0) THEN
      dU = QMIX*Delta_DO
    ELSE
      !write(*,*) ' RESPONSE MIX MIXER_MOD', SCF_ERR
      call dP2(Delta_DO, dU)
      !dU = QMIX*Delta_DO
    ENDIF
    BO = DOrth_old + dU
#endif

    DEALLOCATE(DU, DELTA_DO)

#ifdef PROGRESSON
    CALL BML_DEALLOCATE(AUX_BML)
    CALL BML_DEALLOCATE(DU_BML)
#endif

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
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: X(:,:), YY(:,:), ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: nDelta_DO(:,:), Delta_DS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_1(:,:),H1(:,:), dU(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: MMX(:,:), MMX_I(:,:)

    ALLOCATE(IdentResDM(HDIM,HDIM), Res(NATS), FCOUL_SAVE(3,NATS))
    ALLOCATE(DM_V(HDIM,HDIM,NRANK), UU(HDIM,HDIM,NRANK))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), H_SAVE(HDIM,HDIM))
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
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE, H_SAVE, H_0)
    DEALLOCATE(X, YY, ORTHOH_SAVE, H_1, H1, dU, nDelta_DO, Delta_DS)
  END SUBROUTINE dP2


#ifdef PROGRESSON
  !! progess for dp2, not done yet!
  SUBROUTINE dP2PRG(DELTA_BML, DP_BML)
    IMPLICIT NONE
    TYPE(BML_MATRIX_T), intent(in)  :: DELTA_BML
    TYPE(BML_MATRIX_T), intent(out) :: DP_BML

    !
    INTEGER :: I, J, N, ii, jj
    INTEGER :: KI, MI, K, L
    REAL(LATTEPREC) :: ErrDM, TMP, Nocc, beta, eps, ErrDMn

    REAL(LATTEPREC), ALLOCATABLE :: Res(:), FCOUL_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: DELTAQ_SAVE(:), COULOMBV_SAVE(:)
    REAL(LATTEPREC), ALLOCATABLE :: H_0(:,:), H_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: ORTHOH_SAVE(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: Delta_DS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: H_1(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: MMX(:,:), MMX_I(:,:)

    TYPE(BML_MATRIX_T) :: AUX_BML, SU_BML, DU_BML
    TYPE(BML_MATRIX_T) :: X_BML, Y_BML
    TYPE(BML_MATRIX_T), ALLOCATABLE :: DMV_BML(:), UU_BML(:)

    ALLOCATE(Res(NATS), FCOUL_SAVE(3,NATS))
    ALLOCATE(DELTAQ_SAVE(NATS), COULOMBV_SAVE(NATS))
    ALLOCATE(H_0(HDIM,HDIM), H_SAVE(HDIM,HDIM))
    ALLOCATE(ORTHOH_SAVE(HDIM,HDIM))
    ALLOCATE(Delta_DS(HDIM,HDIM))
    ALLOCATE(H_1(HDIM,HDIM))
    !
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,X_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,Y_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,DP_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,DU_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,SU_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX_BML)

    ALLOCATE(DMV_BML(NRANK),UU_BML(NRANK))
    DO K = 1, NRANK
      CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,DMV_BML(k))
      CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,UU_BML(k))
    ENDDO

    CALL BML_COPY_NEW(DELTA_BML, DU_BML)

    DELTAQ_SAVE   = DELTAQ
    COULOMBV_SAVE = COULOMBV
    FCOUL_SAVE    = FCOUL
    ORTHOH_SAVE   = ORTHOH
    H_SAVE = H
    H_0    = H0
    ErrDM  = 100_LATTEPREC
    eps    = 1.d-7

    ! adpative rank-m, no preconditioning

    do I = 1, HDIM
      do J = 1, HDIM
        SU(I,J) = SMAT(I,J)*DFTB_U(J)
      enddo
    enddo

    CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, SU, SU_BML, ZERO, LT%MDIM)

    K = 0
    ALLOCATE(MMX(1,1),MMX_I(1,1))
    do while ( ErrDM > KERNELTOL .AND. K < NRANK)
      DEALLOCATE(MMX,MMX_I)
      H0 = 0.D0  ! SUCH THAT ADDQDEP ONLY INCLUDES RESPONSE PART
      K = K + 1

      CALL BML_COPY_NEW(DU_BML, DMV_BML(K))

      DO J = 1, K-1
        TMP = BML_SUM_AB(DMV_BML(K),DMV_BML(J),1.D0,1.D-20)
        CALL BML_ADD(DMV_BML(K),DMV_BML(J),1.0D0,-1.D0*TMP,1.D-20) !NUMTHRESH)
      ENDDO

      TMP = SQRT(BML_SUM_SQUARES(DMV_BML(K)))
      IF (TMP<1.d-20) TMP = TMP + 1.d-20
      IF (TMP < 1.D-12) write(6,*) 'TMP is ZERO in dP2!'

      CALL BML_ADD(DMV_BML(K),DMV_BML(K),1.D0/TMP,0.D0,1.D-20) !NUMTHRESH)

      CALL BML_COPY_NEW(ZMAT_BML, AUX_BML)
      CALL BML_MULTIPLY(AUX_BML,DMV_BML(K),Y_BML,1.0D0,0.0D0,NUMTHRESH)
      CALL BML_TRANSPOSE(ZMAT_BML, AUX_BML)
      CALL BML_MULTIPLY(Y_BML,AUX_BML, X_BML,1.0D0,0.0D0,NUMTHRESH)
      CALL BML_MULTIPLY(X_BML,OVER_BML,Y_BML,1.0D0,0.0D0,NUMTHRESH)

      CALL BML_EXPORT_TO_DENSE(Y_BML, Delta_DS)

      CALL BML_TRANSPOSE(Y_BML, AUX_BML)
      CALL BML_MULTIPLY(AUX_BML,SU_BML,X_BML,1.0D0,0.0D0,NUMTHRESH)
      CALL BML_MULTIPLY(SU_BML, X_BML, Y_BML,1.0D0,0.0D0,NUMTHRESH)
      ! H1 stored in X_BML
      CALL BML_ADD(X_BML,Y_BML,-0.25D0,-0.25D0,NUMTHRESH)
      CALL BML_TRANSPOSE(X_BML, AUX_BML)

      ! H_1 stored in X_BML
      CALL BML_ADD(X_BML,AUX_BML,1.0D0,1.0D0,NUMTHRESH)
      CALL BML_EXPORT_TO_DENSE(X_BML, H_1)

      DELTAQ = 0.0
      do I = 1,NATS
        do J = H_INDEX_START(I),H_INDEX_END(I)
          DELTAQ(I) = DELTAQ(I) + 2*Delta_DS(J,J)
        enddo
      enddo

      CALL COULOMBRSPACE
      CALL COULOMBEWALD
      CALL ADDQDEP

      H = H + H_1  ! Total perturbation

      CALL ORTHOMYH
      !CALL ORTHOMYHPRG

      Nocc = BNDFIL*float(HDIM)
      beta = 1.D0/KBT

      ! call canonical response
      call can_resp(ORTHOH,Nocc,beta,EVECS,EVALS,FERMIOCC,CHEMPOT,eps,HDIM)
      !call Canon_DM_PRT(ORTHOH,beta,EVECS,EVALS,CHEMPOT,16,HDIM)

      CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, BO, UU_BML(K), ZERO, LT%MDIM)

      CALL BML_ADD(UU_BML(K),DMV_BML(K),1.0D0,-1.D0,NUMTHRESH)
      CALL BML_COPY_NEW(UU_BML(K), DU_BML)

      !---------------
      !check acccuracy
      !---------------
      ALLOCATE(MMX(K,K),MMX_I(K,K))
      do I = 1,K
        do J = 1,K
          TMP = BML_SUM_AB(UU_BML(I),UU_BML(J),1.D0,1.D-20)
          MMX(I,J) = TMP
        enddo
      enddo

      ! as long as K is small, no need BML here
      call Invert(MMX,MMX_I,K)

      CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX_BML)
      do I = 1,K
        do J = 1,K
          TMP = BML_SUM_AB(DELTA_BML,UU_BML(J),1.D0,1.D-20)
          TMP = TMP * MMX_I(I,J)
          CALL BML_ADD(AUX_BML,UU_BML(I),1.0D0,TMP,NUMTHRESH)
        enddo
      enddo

      ErrDM  = BML_SUM_SQUARES2(AUX_BML,DELTA_BML,1.D0,-1.D0,1.D-20)
      ErrDMn = BML_SUM_SQUARES(AUX_BML)

      IF(ErrDMn<1.d-20) ErrDMn = ErrDMn + 1.d-20
      ErrDM = sqrt(ErrDM)/sqrt(ErrDMn) ! Relative Error
      !write(*,*) ' ErrDM = ', ErrDM  !! Should go down as a function of K
    enddo

    !dP = 0.0_LATTEPREC
    do I = 1,K
      do J = 1,K
        TMP = BML_SUM_AB(UU_BML(J),DELTA_BML,1.D0,1.D-20)
        TMP = - MMX_I(I,J) * TMP

        CALL BML_ADD(DP_BML,DMV_BML(I),1.0D0,TMP,NUMTHRESH)
      enddo
    enddo
    DEALLOCATE(MMX,MMX_I)

    !CALL BML_EXPORT_TO_DENSE(DP_BML, Delta_DS)
    !dP = Delta_DS
    !write(6,*) 'maxval(dP)=', maxval(dabs(dP))

    COULOMBV = COULOMBV_SAVE
    DELTAQ   = DELTAQ_SAVE
    H0 = H_0
    H  = H_SAVE
    FCOUL  = FCOUL_SAVE
    ORTHOH = ORTHOH_SAVE

    DEALLOCATE(Res, FCOUL_SAVE)
    DEALLOCATE(DELTAQ_SAVE, COULOMBV_SAVE, H_SAVE, H_0)
    DEALLOCATE(ORTHOH_SAVE, H_1, Delta_DS)

    CALL BML_DEALLOCATE(X_BML)
    CALL BML_DEALLOCATE(Y_BML)
    !CALL BML_DEALLOCATE(DP_BML)
    CALL BML_DEALLOCATE(DU_BML)
    CALL BML_DEALLOCATE(SU_BML)
    CALL BML_DEALLOCATE(AUX_BML)
    !CALL BML_DEALLOCATE(DELTA_BML)

    DO K = 1, NRANK
      CALL BML_DEALLOCATE(DMV_BML(K))
      CALL BML_DEALLOCATE(UU_BML(K))
    ENDDO
    DEALLOCATE(DMV_BML,UU_BML)

  END SUBROUTINE dP2PRG
#endif
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
      CALL BML_EXPORT_TO_DENSE(DO_BML, DOrth)
      CALL BML_EXPORT_TO_DENSE(DO_BML_OLD, DOrth_old)

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
