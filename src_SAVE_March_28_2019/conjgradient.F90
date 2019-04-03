SUBROUTINE CONJGRADIENT(ITER, DELTAENERGY)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE RELAXCOMMON
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER
  REAL(LATTEPREC) :: DELTAENERGY
  REAL(LATTEPREC) :: MAXSHIFT = 0.1D0, KAPPA = 0.01D0
  REAL(LATTEPREC) :: BETA, BETA1, BETA2, MYSHIFT
  IF (EXISTERROR) RETURN

  ! On the first step, or every so often because
  ! of a loss of conjugacy we do a molecular statics step

  IF ( ITER .EQ. 1 .OR. MOD(ITER, 5) .EQ. 0 ) THEN

     DO I = 1, NATS
        DO J = 1, 3

           !
           ! Limit maximum displacement in case the initial
           ! geometry is not good
           !

           MYSHIFT = KAPPA*FTOT(J,I)

           IF (ABS( MYSHIFT ) .GT. MAXSHIFT) THEN
              MYSHIFT = SIGN(MAXSHIFT, MYSHIFT)
           ENDIF

           CR(J,I) = CR(J,I) + MYSHIFT

        ENDDO
     ENDDO

     OLDF = FTOT
     OLDD = FTOT

  ELSE

     BETA1 = ZERO
     BETA2 = ZERO
     DO I = 1, NATS

        !                                                                    
        ! Use the Polak-Ribero Beta                                          
        !    

        DO J = 1, 3

           BETA1 = BETA1 + FTOT(J,I) * (FTOT(J,I) - OLDF(J,I))
           BETA2 = BETA2 + OLDF(J,I) * OLDF(J,I)

        ENDDO

     ENDDO

     IF ( ABS(BETA2) .LT. 1.0D-6 ) THEN 

        BETA = ZERO

     ELSE

        BETA = BETA1/BETA2

     ENDIF

     DO I = 1, NATS

        DO J = 1, 3

           D1(J,I) = FTOT(J,I) + BETA*OLDD(J,I)

           MYSHIFT = KAPPA*D1(J,I)

           IF ( ABS( MYSHIFT ) .GT. MAXSHIFT ) THEN
              MYSHIFT = SIGN( MAXSHIFT, MYSHIFT )
           ENDIF

           CR(J,I) = CR(J,I) + MYSHIFT

        ENDDO
     ENDDO

     OLDF = FTOT
     OLDD = D1

  ENDIF

  RETURN

END SUBROUTINE CONJGRADIENT
