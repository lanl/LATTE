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

SUBROUTINE XBODM(ITER)  !! ANDERS CHANGE NEW SUBROUTINE FOR DM PROPAGATION FOR DFTB+U WITH XLBOMD

  !
  ! This here subroutine implements Niklasson's eXtended
  ! Born-Oppenheimer, time-reversible, and most excellent MD thing.
  !
  ! Niklasson, PRL, vol. 100, p 123004, (2008)
  !

  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SETUPARRAY
  USE MYPRECISION
  USE NONOARRAY
  USE DMARRAY

!             DELTAQ(I) = -ATOCC(ELEMPOINTER(I))

  IMPLICIT NONE
  REAL(LATTEPREC) :: D(HDIM,HDIM)

  INTEGER :: I, J, K, ITER

  IF (EXISTERROR) RETURN

  !
  ! If we have no damping:
  !
  ! Our time reversible guess for the next set of diagonal H-matrix
  ! elements depend on those from the last guess, the guess before that, and
  ! the H-matrix elements we calculated subject to contraints on Delta_q
  ! from the last iteration.
  !
  ! P(t) - last guess
  ! P(t - dt) - guess before that
  ! D(t) - the elements calculated from the last run
  !
  ! P(t + dt) = 2P(t) - P(t - dt) + 2[D(t) - P(t)]
  !
  ! With damping, it's the same general idea but we use guesses from even
  ! earlier iterations too. Preliminary testing shows better energy
  ! conservation with high 'K'.
  !

  ! Required for the Fast-QMMD scheme

  IF (QITER .EQ. 0) THEN
     KAPPA_SCALE = MDMIX
  ELSE
     KAPPA_SCALE = ONE
  ENDIF

  IF (ITER .EQ. 1) THEN

     IF (ELECTRO .EQ. 1) THEN

       ! DO Orthogonal DM representation
       PO_0 = DOrth
       PO_1 = DOrth
       PO_2 = DOrth
       PO_3 = DOrth
       PO_4 = DOrth
       PO_5 = DOrth
       PO_6 = DOrth
       PO_7 = DOrth
       POrth = DOrth

        call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,POrth,HDIM,ZERO,NONOTMP,HDIM)
        call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,NONOTMP,HDIM,XMAT,HDIM,ZERO,PNO,HDIM)
        D = PNO
        do I = 1, NATS
          DELTAQDM(I) = -ATOCC(ELEMPOINTER(I))
          do K = H_INDEX_START(I), H_INDEX_END(I)
              DELTAQDM(I) = DELTAQDM(I) + dot_product(PNO(:,K),SMAT(:,K))
          enddo
        enddo

     ELSEIF (ELECTRO .EQ. 0) THEN

        WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
        STOP

     ENDIF

  ELSEIF (ITER .GT. 1) THEN

     IF (XBODISON .EQ. 0) THEN

        IF (ELECTRO .EQ. 0) THEN

          WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
          STOP

        ELSEIF (ELECTRO .EQ. 1) THEN

              POrth = TWO*PO_0 - PO_1 + TWO*(DOrth-POrth) 
              PO_1 = PO_0; PO_0 = POrth;
                          
        ENDIF


     ELSEIF (XBODISON .EQ. 1) THEN

        IF (ELECTRO .EQ. 0) THEN

          WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
          STOP

        ELSEIF (ELECTRO .EQ. 1) THEN

            IF (ITER <= 2) THEN  ! NEEDS TO BE INITIATE THE d2PO CALCULATIONS
               POrth = 2.D0*PO_0 - PO_1 + KAPPA_SCALE*KAPPA_XBO*(DOrth-PO_0) &  ! OLD VERSION Kernel = c*I
                   + ALPHA_XBO*(CNK(1)*PO_0+CNK(2)*PO_1+CNK(3)*PO_2+CNK(4)*PO_3+CNK(5)*PO_4 &
                   + CNK(6)*PO_5+CNK(7)*PO_6)
            ELSE
                POrth = 2.D0*PO_0 - PO_1 + KAPPA_XBO*d2PO &  !! NEW DM PROPAGATION WITH RANK-1 UPDATES
                   + ALPHA_XBO*(CNK(1)*PO_0+CNK(2)*PO_1+CNK(3)*PO_2+CNK(4)*PO_3+CNK(5)*PO_4 &
                   + CNK(6)*PO_5+CNK(7)*PO_6)
            ENDIF
            PO_6 = PO_5; PO_5 = PO_4; PO_4 = PO_3; PO_3 = PO_2; 
            PO_2 = PO_1; PO_1 = PO_0; PO_0 = POrth;

            IF (FULLQCONV==1) POrth = DOrth

            call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,POrth,HDIM,ZERO,NONOTMP,HDIM)
            call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,NONOTMP,HDIM,XMAT,HDIM,ZERO,PNO,HDIM)

           do I = 1, NATS
             DELTAQDM(I) = 0.D0
             do K = H_INDEX_START(I), H_INDEX_END(I)
                 DELTAQDM(I) = DELTAQDM(I) + dot_product(PNO(:,K),SMAT(:,K))
             enddo
             DELTAQDM(I) = DELTAQDM(I) - ATOCC(ELEMPOINTER(I))
           enddo

           !!! OBS DELTAQDM IS THE NEW UPDATED PARTIAL OCCUPATION, NOT DELTAQ THAT NEEDS TO BE UPDATED SEPARATELY

        ENDIF
     ENDIF
  ENDIF


  RETURN

END SUBROUTINE XBODM



SUBROUTINE XBODMPRG(ITER)  !! 
  ! This here subroutine implements Niklasson's eXtended
  ! Born-Oppenheimer, time-reversible, and most excellent MD thing.
  !
  ! Niklasson, PRL, vol. 100, p 123004, (2008)
  !
#ifdef PROGRESSON
  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SETUPARRAY
  USE MYPRECISION
  USE NONOARRAY
  USE DMARRAY
  USE SPARSEARRAY
  USE BML
  USE LATTEPARSER
  USE GENXPROGRESS
  USE NONOARRAYPROGRESS
!             DELTAQ(I) = -ATOCC(ELEMPOINTER(I))
  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = LATTEPREC
  TYPE(BML_MATRIX_T) :: AUX1_BML
  TYPE(BML_MATRIX_T) :: AUX2_BML

  INTEGER :: I, J, K, ITER

  IF (EXISTERROR) RETURN

  !
  ! If we have no damping:
  !
  ! Our time reversible guess for the next set of diagonal H-matrix
  ! elements depend on those from the last guess, the guess before that, and
  ! the H-matrix elements we calculated subject to contraints on Delta_q
  ! from the last iteration.
  !
  ! P(t) - last guess
  ! P(t - dt) - guess before that
  ! D(t) - the elements calculated from the last run
  !
  ! P(t + dt) = 2P(t) - P(t - dt) + 2[D(t) - P(t)]
  !
  ! With damping, it's the same general idea but we use guesses from even
  ! earlier iterations too. Preliminary testing shows better energy
  ! conservation with high 'K'.
  !

  ! Required for the Fast-QMMD scheme

  IF (QITER .EQ. 0) THEN
     KAPPA_SCALE = MDMIX
  ELSE
     KAPPA_SCALE = ONE
  ENDIF

  CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX1_BML)
  CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX2_BML)
  !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, DOrth, DO_BML, ZERO, LT%MDIM)

  WRITE(6,*) 'test-xbodm: current iter', ITER, XBODISON
  IF (ITER .EQ. 1) THEN

     IF (ELECTRO .EQ. 1) THEN

       ! DO Orthogonal DM representation
       CALL BML_COPY_NEW(DO_BML,PO0_BML)
       CALL BML_COPY_NEW(DO_BML,PO1_BML)
       CALL BML_COPY_NEW(DO_BML,PO2_BML)
       CALL BML_COPY_NEW(DO_BML,PO3_BML)
       CALL BML_COPY_NEW(DO_BML,PO4_BML)
       CALL BML_COPY_NEW(DO_BML,PO5_BML)
       CALL BML_COPY_NEW(DO_BML,PO6_BML)
       CALL BML_COPY_NEW(DO_BML,PO_BML)

       
       CALL BML_MULTIPLY(ZMAT_BML,PO_BML,AUX1_BML,1.0_dp,0.0_dp,NUMTHRESH)
       CALL BML_TRANSPOSE(ZMAT_BML, AUX2_BML)
       CALL BML_MULTIPLY(AUX1_BML,AUX2_BML,PNO_BML,1.0_dp,0.0_dp,NUMTHRESH)
       CALL BML_EXPORT_TO_DENSE(PNO_BML, PNO)

      ! call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,POrth,HDIM,ZERO,NONOTMP,HDIM)
      ! call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,NONOTMP,HDIM,XMAT,HDIM,ZERO,PNO,HDIM)

       do I = 1, NATS
         DELTAQDM(I) = -ATOCC(ELEMPOINTER(I))
         do K = H_INDEX_START(I), H_INDEX_END(I)
             DELTAQDM(I) = DELTAQDM(I) + dot_product(PNO(:,K),SMAT(:,K))
         enddo
       enddo

     ELSEIF (ELECTRO .EQ. 0) THEN

        WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
        STOP

     ENDIF

  ELSEIF (ITER .GT. 1) THEN

     IF (XBODISON .EQ. 0) THEN

        IF (ELECTRO .EQ. 0) THEN

          WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
          STOP

        ELSEIF (ELECTRO .EQ. 1) THEN
              
              POrth = TWO*PO_0 - PO_1 + TWO*(DOrth-POrth) 
              PO_1 = PO_0; PO_0 = POrth;
                          
              CALL BML_ADD(PO_BML,DO_BML,-1.0_DP,1.0_DP,NUMTHRESH)
              CALL BML_ADD(PO_BML,PO1_BML,2.0_DP,-1.0_DP,NUMTHRESH)
              CALL BML_ADD(PO_BML,PO0_BML,1.0_DP, 2.0_DP,NUMTHRESH)

              CALL BML_COPY_NEW(PO0_BML,PO1_BML)
              CALL BML_COPY_NEW(PO_BML,PO0_BML)
        ENDIF


     ELSEIF (XBODISON .EQ. 1) THEN

        IF (ELECTRO .EQ. 0) THEN

          WRITE(*,*) ' ELECTRO = 0 Not implemented, stop '
          STOP

        ELSEIF (ELECTRO .EQ. 1) THEN

            IF (ITER <= 2) THEN  ! NEEDS TO BE INITIATE THE d2PO CALCULATIONS

               !POrth = 2.D0*PO_0 - PO_1 + KAPPA_SCALE*KAPPA_XBO*(DOrth-PO_0) &  ! OLD VERSION Kernel = c*I
               !    + ALPHA_XBO*(CNK(1)*PO_0+CNK(2)*PO_1+CNK(3)*PO_2+CNK(4)*PO_3+CNK(5)*PO_4 &
               !    + CNK(6)*PO_5+CNK(7)*PO_6)

               CALL BML_COPY_NEW(PO6_BML,AUX1_BML)

               CALL BML_ADD(AUX1_BML,PO5_BML,CNK(7),CNK(6),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO4_BML,1.0_DP,CNK(5),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO3_BML,1.0_DP,CNK(4),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO2_BML,1.0_DP,CNK(3),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO1_BML,1.0_DP,CNK(2),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO0_BML,1.0_DP,CNK(1),NUMTHRESH)

               CALL BML_ADD(AUX1_BML,PO0_BML,ALPHA_XBO,2.0_DP,NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO1_BML,1.0_DP,-1.0_DP,NUMTHRESH)

               CALL BML_COPY_NEW(DO_BML,PO_BML)
               CALL BML_ADD(PO_BML,PO0_BML,1.0_DP,-1.0_DP,NUMTHRESH)
               CALL BML_ADD(PO_BML,AUX1_BML,KAPPA_SCALE*KAPPA_XBO,1.0_DP,NUMTHRESH)
               
            ELSE
               !POrth = 2.D0*PO_0 - PO_1 + KAPPA_XBO*d2PO &  !! NEW DM PROPAGATION WITH RANK-1 UPDATES
               !    + ALPHA_XBO*(CNK(1)*PO_0+CNK(2)*PO_1+CNK(3)*PO_2+CNK(4)*PO_3+CNK(5)*PO_4 &
               !    + CNK(6)*PO_5+CNK(7)*PO_6)

               CALL BML_COPY_NEW(PO6_BML,AUX1_BML)
               CALL BML_ADD(AUX1_BML,PO5_BML,CNK(7),CNK(6),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO4_BML,1.0_DP,CNK(5),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO3_BML,1.0_DP,CNK(4),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO2_BML,1.0_DP,CNK(3),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO1_BML,1.0_DP,CNK(2),NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO0_BML,1.0_DP,CNK(1),NUMTHRESH)

               CALL BML_ADD(AUX1_BML,PO0_BML,ALPHA_XBO,2.0_DP,NUMTHRESH)
               CALL BML_ADD(AUX1_BML,PO1_BML,1.0_DP,-1.0_DP,NUMTHRESH)

               CALL BML_COPY_NEW(D2P_BML, PO_BML)
               !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, D2PO, PO_BML, ZERO, LT%MDIM)
               CALL BML_ADD(PO_BML,AUX1_BML,KAPPA_XBO,1.0_DP,NUMTHRESH)
            ENDIF
            !PO_6 = PO_5; PO_5 = PO_4; PO_4 = PO_3; PO_3 = PO_2; 
            !PO_2 = PO_1; PO_1 = PO_0; PO_0 = POrth;
            CALL BML_COPY_NEW(PO5_BML,PO6_BML)
            CALL BML_COPY_NEW(PO4_BML,PO5_BML)
            CALL BML_COPY_NEW(PO3_BML,PO4_BML)
            CALL BML_COPY_NEW(PO2_BML,PO3_BML)
            CALL BML_COPY_NEW(PO1_BML,PO2_BML)
            CALL BML_COPY_NEW(PO0_BML,PO1_BML)
            CALL BML_COPY_NEW(PO_BML,PO0_BML)

            !IF (FULLQCONV==1) POrth = DOrth
            IF (FULLQCONV==1) CALL BML_COPY_NEW(DO_BML,PO_BML)


            CALL BML_MULTIPLY(ZMAT_BML,PO_BML,AUX1_BML,1.0_dp,0.0_dp,NUMTHRESH)
            CALL BML_TRANSPOSE(ZMAT_BML, AUX2_BML)
            CALL BML_MULTIPLY(AUX1_BML,AUX2_BML,PNO_BML,1.0_dp,0.0_dp,NUMTHRESH)
            CALL BML_EXPORT_TO_DENSE(PNO_BML, PNO)

            !call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,XMAT,HDIM,POrth,HDIM,ZERO,NONOTMP,HDIM)
            !call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,NONOTMP,HDIM,XMAT,HDIM,ZERO,PNO,HDIM)

           do I = 1, NATS
             DELTAQDM(I) = 0.D0
             do K = H_INDEX_START(I), H_INDEX_END(I)
                 DELTAQDM(I) = DELTAQDM(I) + dot_product(PNO(:,K),SMAT(:,K))
             enddo
             DELTAQDM(I) = DELTAQDM(I) - ATOCC(ELEMPOINTER(I))
           enddo

           !!! OBS DELTAQDM IS THE NEW UPDATED PARTIAL OCCUPATION, NOT DELTAQ THAT NEEDS TO BE UPDATED SEPARATELY

        ENDIF
     ENDIF
  ENDIF

  CALL BML_DEALLOCATE(AUX1_BML)
  CALL BML_DEALLOCATE(AUX2_BML)

  RETURN

#endif
END SUBROUTINE XBODMPRG
