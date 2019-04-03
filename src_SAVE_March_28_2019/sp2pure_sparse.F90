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

SUBROUTINE SP2PURE_SPARSE

  !
  ! This subroutine implements Niklasson's SP2 density matrix purification
  ! algorithm.
  !

#ifdef DBCSR_ON
  !**************************************************************
  !use these dbcsr mod files
  USE DBCSR_VAR_MOD
  USE dbcsr_types
  USE dbcsr_methods
  USE dbcsr_error_handling
  USE array_types,                     ONLY: array_data,&
       array_i1d_obj,&
       array_new,&
       array_nullify,&
       array_release,&
       array_size
  USE dbcsr_io
  USE dbcsr_operations
  USE dbcsr_ptr_util
  USE dbcsr_transformations
  USE dbcsr_util
  USE dbcsr_work_operations
  USE dbcsr_message_passing

  USE dbcsr_block_access
  USE dbcsr_iterator_operations,       ONLY: dbcsr_iterator_blocks_left,&
       dbcsr_iterator_next_block,&
       dbcsr_iterator_start,&
       dbcsr_iterator_stop

  USE dbcsr_dist_operations,          ! ONLY: create_bl_distribution,&                                                               !        dbcsr_get_stored_coordinates
  !*****************************************************************
#endif

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPARSEARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, ITERZ, a, b
  INTEGER :: BREAKLOOP
  INTEGER :: IP, JP, KP, NNZX2, PREVNNZX2
  INTEGER :: NEWROW
  INTEGER :: COUNT, KEEP, PREVKEEP
  INTEGER :: NEWXB

  !
  ! Naming: CX = ColIndX  RX = RowIndX   VX = ValX
  !

  INTEGER, ALLOCATABLE :: CX(:), CXTMP(:)
  REAL(LATTEPREC) :: TRX, OCC
  REAL(LATTEPREC) :: TRX2, TR2XX2, TRXOLD
  REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0E-14
  REAL(LATTEPREC), ALLOCATABLE :: VX(:), VXTMP(:)
  REAL(LATTEPREC) :: LIMDIFF, ABSHAM
  REAL(LATTEPREC) :: CIRCLE(2), RADIUS, GERSHFACT
  REAL(LATTEPREC), ALLOCATABLE :: BLOCK_2D(:,:)
  REAL(LATTEPREC) :: FROB

  !
  ! We're also using Niklasson's scheme to determine convergence
  !

  REAL(LATTEPREC) :: IDEMPERR = ZERO, IDEMPERR1 = ZERO, IDEMPERR2 = ZERO
  IF (EXISTERROR) RETURN

  OCC = BNDFIL*FLOAT(HDIM)

  TRX = ZERO

  ! build the sparse H matrix and compute the bounds using
  ! Gershgorin circles

  CALL GERSHGORIN

#ifdef DBCSR_ON

  IF (SPINON .EQ. 1) THEN
     CALL ERRORS("sp2pure_sparse","Sparse SP2 with DBSCR for &
          & spin-polarized systems hasn't been implemented just yet")
  ENDIF


  !
  ! We have to put things into an array which has been padded
  ! to handle the block size
  !

  BO_PADDED = ZERO

  IF (BASISTYPE .EQ. "ORTHO") THEN

     DO I = 1, HDIM
        DO J = 1, HDIM
           BO_PADDED(J,I) = -H(J,I)/MAXMINUSMIN
        ENDDO
     ENDDO

  ELSE

     DO I = 1, HDIM
        DO J = 1, HDIM
           BO_PADDED(J,I) = -ORTHOH(J,I)/MAXMINUSMIN
        ENDDO
     ENDDO

  ENDIF

  GERSHFACT =  MAXEVAL/MAXMINUSMIN


  !important for quick convergence

  TRX = ZERO

  DO I = 1, HDIM

     BO_PADDED(I,I) = GERSHFACT + BO_PADDED(I,I)
     TRX = TRX + BO_PADDED(I,I)

  ENDDO

  !initiallizes matrix a

  CALL dbcsr_init(matrix_a)
  CALL dbcsr_init(matrix_b)

  !create the dbcsr matrix, a double precision non symmetric matrix
  !with nblkrows_total x nblkcols_total blocks and
  !sizes "sum(row_blk_sizes)" x "sum(col_blk_sizes)", distributed as
  !specified by the dist_a object

  CALL dbcsr_create (matrix_a,"Matrix A ",&
       dist_a, dbcsr_type_no_symmetry,&
       row_blk_sizes,&
       col_blk_sizes,&
       data_type=dbcsr_type_real_8,&
       error=error)

  CALL dbcsr_create (matrix_b,"Matrix B ",&
       dist_a, dbcsr_type_no_symmetry,&
       row_blk_sizes,&
       col_blk_sizes,&
       data_type=dbcsr_type_real_8,&
       error=error)

  !puts values into matrix_a

  ALLOCATE(BLOCK_2D(BLKSZ,BLKSZ))

  DO j=1, dbcsr_nblkcols_total(matrix_a)
     DO i=1, dbcsr_nblkrows_total(matrix_a)

	FROB = ZERO
        DO a = 1, BLKSZ
           DO b = 1, BLKSZ

              BLOCK_2D(b,a)=BO_PADDED(BLKSZ*(i-1)+b,BLKSZ*(j-1)+a)
              FROB = FROB + BLOCK_2D(b,a)*BLOCK_2D(b,a)

           ENDDO
        ENDDO

        tr = .FALSE.

        IF (FROB .GT. 1.0D-12) THEN


           CALL dbcsr_get_stored_coordinates (matrix_a, i, j, tr, proc_holds_blk)

           IF(proc_holds_blk .EQ. dbcsr_mp_mynode(dbcsr_distribution_mp (dbcsr_distribution (matrix_a)))) THEN

              !call to put a block into the matrix

              CALL dbcsr_put_block(matrix_a, i, j, BLOCK_2D)

           ENDIF
        ENDIF
     ENDDO
  ENDDO


  !refinallizes matrix a

  CALL dbcsr_finalize(matrix_a, error=error)
  CALL dbcsr_finalize(matrix_b, error=error)

  ! Now throw away all the BLKSZ*BLKSZ blocks that have nothing in them

  ! OCC = dbcsr_get_occupation(matrix_a)
  ! PRINT*, "OCC 1 = ", OCC

  ! CALL dbcsr_filter(matrix_a, eps=1.0d-12, error=error)

  ! OCC = dbcsr_get_occupation(matrix_a)
  ! PRINT*, "OCC 2 = ", OCC

  ITERZ = 0

  BREAKLOOP = 0

  DO WHILE ( BREAKLOOP .EQ. 0 .AND. ITERZ .LT. 100 )

     ITERZ = ITERZ + 1

     CALL dbcsr_copy(matrix_b, matrix_a, error=error)

     IF (THRESHOLDON .EQ. 1) THEN

	! Since we threshold based on the Frobenius norm, we need to multiply
 	! the numerical threshold by the block size in order to maintain a
	! connection with our element-by-element algorithm
	! X_f = sqrt(sum a_ij^2) -> sqrt(blksz ^2 * eps^2) = blksz*eps

        CALL dbcsr_multiply('n','n', MINUSONE, matrix_a, matrix_a, ONE,&
             matrix_b, filter_eps=FLOAT(BLKSZ)*NUMTHRESH, error=error)

     ELSE

        CALL dbcsr_multiply('n','n', MINUSONE, matrix_a, matrix_a, ONE,&
             matrix_b, error=error)

     ENDIF

     TRX2 = ZERO

     CALL dbcsr_trace(matrix_b, TRX2, error=error)

     LIMDIFF = ABS(TRX - TRX2 - OCC) - ABS(TRX + TRX2 - OCC)

     IF ( LIMDIFF .GE. IDEMTOL ) THEN

        !X <- X + X_temp
        CALL dbcsr_add(matrix_a, matrix_b, ONE, ONE, error=error)

        TRX = TRX + TRX2

     ELSEIF ( LIMDIFF .LT. -IDEMTOL ) THEN

        !X <- X - X_temp
        CALL dbcsr_add(matrix_a, matrix_b, ONE, MINUSONE, error=error)

        TRX = TRX - TRX2

     ENDIF

     ! Filtering, if desired

     !     IF (THRESHOLDON .EQ. 1) CALL dbcsr_filter(matrix_a, &
     !          eps=NUMTHRESH, error=error)

     IDEMPERR2 = IDEMPERR1
     IDEMPERR1 = IDEMPERR
     IDEMPERR = ABS(TRX2)


     !	WRITE(*,10) ITERZ, IDEMPERR, IDEMPERR2 - IDEMPERR
10   FORMAT(I4, 2G30.18)
     IF (SP2CONV .EQ. "REL" .AND. ITERZ .GE. MINSP2ITER &
          .AND. (IDEMPERR2 .LE. IDEMPERR .OR. &
          IDEMPERR .LT. IDEMTOL)) BREAKLOOP = 1

     !	IF (ITERZ .EQ. 30) BREAKLOOP=1
     IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LT. IDEMTOL) BREAKLOOP = 1

  ENDDO

  IF (ITERZ .EQ. 100) THEN
     CALL PANIC
     CALL ERRORS("sp2pure_sparse","SP2 purification is not converging: STOP!")
  ENDIF


  !shares data with all procs
  CALL dbcsr_replicate_all(matrix_a, error)

  !  BO = ZERO


  DO j=1, dbcsr_nblkcols_total(matrix_a)
     DO i=1, dbcsr_nblkrows_total(matrix_a)

        tr = .FALSE.

        CALL dbcsr_get_stored_coordinates (matrix_a, i, j, tr, proc_holds_blk)

        !call to put a block into the matrix

        CALL dbcsr_get_block(matrix_a, i, j, BLOCK_2D, tr, found)

	DO a = 1, BLKSZ
           DO b = 1, BLKSZ

              IF(found) THEN
                 BO_PADDED(BLKSZ*(i-1)+b,BLKSZ*(j-1)+a) = TWO*BLOCK_2D(b,a)
              ENDIF

           ENDDO
	ENDDO


     ENDDO
  ENDDO
  DEALLOCATE(BLOCK_2D)


  CALL dbcsr_release(matrix_a)
  CALL dbcsr_release(matrix_b)

  DO I = 1, HDIM
     DO J = 1, HDIM
        BO(J,I) = BO_PADDED(J,I)
     ENDDO
  ENDDO


#elif defined(DBCSR_OFF)


  ALLOCATE( VX(NNZ), CX(NNZ) )

  COUNT = 0

  ! Put the normalized H into compressed sparse row format

  IF (BASISTYPE .EQ. "ORTHO") THEN

     DO I = 1, HDIM

        NEWROW = 0

        DO J = 1, HDIM

           IF ( ABS( H(J,I) ) .GT. HTHRESH ) THEN

              COUNT = COUNT + 1

              IF (I .EQ. J) THEN

                 VX( COUNT ) = (MAXEVAL - H(J,I))/MAXMINUSMIN

              ELSE

                 VX( COUNT ) = -H(J,I)/MAXMINUSMIN

              ENDIF

              CX( COUNT ) = J

              IF ( NEWROW .EQ. 0 ) THEN
                 RX( I ) = COUNT
                 NEWROW = 1
              ENDIF

           ENDIF

        ENDDO
     ENDDO

  ELSE

     DO I = 1, HDIM

        NEWROW = 0

        DO J = 1, HDIM

           IF ( ABS( ORTHOH(J,I) ) .GT. HTHRESH ) THEN

              COUNT = COUNT + 1

              IF (I .EQ. J) THEN

                 VX( COUNT ) = (MAXEVAL - ORTHOH(J,I))/MAXMINUSMIN

              ELSE

                 VX( COUNT ) = -ORTHOH(J,I)/MAXMINUSMIN

              ENDIF

              CX( COUNT ) = J

              IF ( NEWROW .EQ. 0 ) THEN
                 RX( I ) = COUNT
                 NEWROW = 1
              ENDIF

           ENDIF

        ENDDO
     ENDDO

  ENDIF

  RX(HDIM + 1) = COUNT + 1

  ITERZ = 0
  BREAKLOOP = 0

  PREVNNZX2 = 0
  PREVKEEP = 0

  DO WHILE ( BREAKLOOP .EQ. 0 .AND. ITERZ .LT. 100 )

     ITERZ = ITERZ + 1

     !
     ! Sparse BO * BO
     !

     ! First we need to determine the number of non-zeros in X*X
     ! so we can allocate storage


     IF (ITERZ .LE. FILLINSTOP) THEN

        IP = 0

        DO I = 1, HDIM

           XB = 0

           DO JP = RX( I ), RX( I + 1 ) - 1

              J = CX( JP )

              DO KP = RX( J ), RX( J + 1 ) - 1

                 K = CX( KP )

                 IF ( XB( K ) .NE. I ) THEN

                    IP = IP + 1
                    XB( K ) = I

                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        NNZX2 = IP

     ENDIF

     ! Allocate a bit more storage just to be safe

     IF (ITERZ .EQ. FILLINSTOP) NNZX2 = INT(1.1*NNZX2)

     !     print*, ITERZ, NNZX2, PREVNNZX2, FLOAT(NNZX2)/FLOAT((HDIM*HDIM))

     IF (ITERZ .EQ. 1) THEN

        ALLOCATE( VXTMP( NNZX2 ), CXTMP( NNZX2 ))

     ELSEIF (ITERZ .GT. 1 .AND. NNZX2 .GT. PREVNNZX2) THEN

        ! We will not enter for ITERZ > FILLINSTOP because
        ! of the condition PREVNNZX2 = NNZX2

        DEALLOCATE( VXTMP, CXTMP )
        ALLOCATE( VXTMP( NNZX2 ), CXTMP( NNZX2 ) )

        PREVNNZX2 = NNZX2

     ENDIF

     IP = 1
     XB = 0

     TR2XX2 = ZERO
     TRX2 = ZERO

     DO I = 1, HDIM

        RXTMP(I) = IP

        DO JP = RX( I ), RX( I + 1 ) - 1

           J = CX( JP )

           ! Computing these traces in this loop rather than separately

           IF (J .EQ. I) TR2XX2 = TR2XX2 + VX(JP)

           DO KP = RX( J ), RX( J + 1 ) - 1

              K = CX( KP )

              IF ( XB( K ) .NE. I ) THEN

                 CXTMP( IP ) = K
                 IP = IP + 1
                 XB( K ) = I
                 WORK( K ) = VX( JP ) * VX( KP )

              ELSE

                 WORK( K ) = WORK( K ) + VX( JP ) * VX( KP )

              ENDIF

           ENDDO

        ENDDO

        DO J = RXTMP( I ), IP - 1

           VXTMP( J ) = WORK( CXTMP( J ) )

           ! same here

           IF (CXTMP(J) .EQ. I) TRX2 = TRX2 + VXTMP( J )

        ENDDO

     ENDDO

     RXTMP( HDIM + 1 ) = IP

     TR2XX2 = TWO*TR2XX2 - TRX2

     TRXOLD = TRX

     LIMDIFF = ABS(TRX2 - OCC) - ABS(TR2XX2 - OCC)

     IF ( LIMDIFF .LT. -IDEMTOL ) THEN

        !
        ! X <= X^2
        !

        TRX = TRX2

     ELSEIF ( LIMDIFF .GT. IDEMTOL ) THEN

        TRX = TR2XX2

        !
        ! X <= 2X - X^2
        !

        DO I = 1, HDIM

           DO IP = RXTMP( I ), RXTMP( I + 1 ) - 1

              WORK( CXTMP( IP ) ) = -VXTMP( IP )

           ENDDO

           DO IP = RX( I ), RX( I + 1 ) - 1

              WORK( CX( IP ) ) = WORK( CX( IP ) ) + TWO*VX( IP )

           ENDDO

           DO IP = RXTMP( I ), RXTMP( I + 1 ) - 1

              VXTMP( IP ) = WORK( CXTMP( IP ) )

           ENDDO

        ENDDO

     ENDIF

     IF (ABS(LIMDIFF) .GT. IDEMTOL) THEN

        IF (THRESHOLDON .EQ. 1) THEN


           IF (ITERZ .LT. FILLINSTOP) THEN

              KEEP = 0
              DO I = 1, NNZX2
                 IF ( ABS( VXTMP( I ) ) .GT. NUMTHRESH ) KEEP = KEEP + 1
              ENDDO

              IF ( KEEP .GT. PREVKEEP ) THEN
                 DEALLOCATE( VX, CX )
                 ALLOCATE( VX( KEEP ), CX( KEEP ) )
                 PREVKEEP = KEEP
              ENDIF

           ELSEIF (ITERZ .EQ. FILLINSTOP) THEN

              DEALLOCATE( VX, CX )
              ALLOCATE( VX( NNZX2 ), CX( NNZX2 ) )

           ENDIF

           COUNT = 0
           DO I = 1, HDIM
              NEWROW = 0
              DO J = RXTMP( I ), RXTMP( I + 1 ) - 1

                 ! Throwing away small elements while copying

                 IF ( ABS( VXTMP( J ) ) .GT. NUMTHRESH ) THEN

                    COUNT = COUNT + 1
                    VX( COUNT ) = VXTMP( J )
                    CX( COUNT ) = CXTMP( J )

                    IF ( NEWROW .EQ. 0 ) THEN
                       RX( I ) = COUNT
                       NEWROW = 1
                    ENDIF

                 ENDIF
              ENDDO
           ENDDO

           RX( HDIM + 1 ) = COUNT + 1

        ELSE

           IF (ITERZ .LE. FILLINSTOP .AND. RX(HDIM + 1) .NE. NNZX2 ) THEN

              DEALLOCATE( VX, CX )
              ALLOCATE( VX(NNZX2), CX(NNZX2) )

           ENDIF

           VX = VXTMP
           CX = CXTMP
           RX = RXTMP

        ENDIF

     ENDIF

     IDEMPERR2 = IDEMPERR1
     IDEMPERR1 = IDEMPERR
     IDEMPERR = ABS(TRX - TRXOLD)

     !     PRINT*, IDEMPERR, ABS(TRX2 - OCC) - ABS(TR2XX2 - OCC)
     !     PRINT*, ITERZ, IDEMPERR0 - IDEMPERR2, ABS(TRX - OCC)
     !PRINT*, ITERZ, ABS(TRX - OCC)

     IF (SP2CONV .EQ. "REL" .AND.  ITERZ .GE. MINSP2ITER &
          .AND. (IDEMPERR .GE. IDEMPERR2) ) BREAKLOOP = 1

     IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LE. IDEMTOL) BREAKLOOP = 1

  ENDDO

  IF (ITERZ .EQ. 100) THEN
     CALL PANIC
     CALL ERRORS("sp2pure_sparse","Sparse SP2 purification is not converging: STOP!")
  ENDIF

  BO = ZERO

  DO I = 1, HDIM
     DO J = RX( I ), RX( I + 1 ) - 1

        BO( I, CX(J) ) = TWO*VX( J )

     ENDDO
  ENDDO

  DEALLOCATE(VX, VXTMP, CX, CXTMP)

#endif

  RETURN

END SUBROUTINE SP2PURE_SPARSE
