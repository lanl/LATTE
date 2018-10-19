PROGRAM SPLINE

  ! This is taken from Numerical recipes
  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER, PARAMETER :: N = 499
  REAL :: X(N), Y(N), Y2(N), P, QN, SIG, UN, U(N)
  REAL :: R, NEWY,GRAD

  OPEN(UNIT=11, STATUS="OLD", FILE="inputpoints.dat")
  OPEN(UNIT=12, STATUS="UNKNOWN", FILE="checknumrep.dat")

  DO I = 1, N
     READ(11,*) X(I), Y(I)
  ENDDO

  Y2(1) = 0.0
  U(1) = 0.0

  DO I = 2, N-1
     SIG = (X(I) - X(I-1))/(X(I+1) - X(I-1))
     P = SIG*Y2(I-1) + 2.0
     Y2(I) = (SIG - 1.0)/P
     U(I) = (6.0*((Y(I+1) - Y(I))/(X(I+1) - X(I)) - (Y(I) - Y(I-1)) &
          /(X(I) - X(I-1)))/(X(I+1)-X(I-1)) - SIG*U(I-1))/P
  ENDDO

  QN = 0.0
  UN = 0.0

  Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1) + 1.0)

  DO K = N-1, 1, -1
     Y2(K) = Y2(K)*Y2(K+1) + U(K)
  ENDDO

  DO I = 1, 1000

     R = 0.5 + 0.5*REAL(I-1)/1000.0


     CALL SPLINT(X, Y, Y2, R, NEWY, GRAD)

     WRITE(12,*) R, NEWY, GRAD
  ENDDO

END PROGRAM SPLINE

SUBROUTINE SPLINT(X, Y, Y2, R, NEWY,GRAD)

  IMPLICIT NONE

  INTEGER :: K, KHI, KLO
  INTEGER, PARAMETER :: N = 499
  REAL :: X(N), Y(N), Y2(N), R, NEWY, A, B, H, GRAD

  KLO = 1
  KHI = N

  DO WHILE(KHI - KLO .GT. 1)
     K = (KHI + KLO)/2
     IF (X(K) .GT. R) THEN
        KHI = K
     ELSE
        KLO = K
     ENDIF
  ENDDO

  H = X(KHI) - X(KLO)

  A = (X(KHI) - R)/H
  B = (R - X(KLO))/H

  NEWY = A*Y(KLO) + B*Y(KHI) + &
       ((A*A*A - A)*Y2(KLO) + (B*B*B - B)*Y2(KHI))*(H*H/6.0)
  
  GRAD = -Y(KLO)/H + Y(KHI)/H + &
       ((3.0*A*A*(-1.0/H) + 1.0/H)*Y2(KLO) + &
       (3.0*B*B*(1.0/H) - 1.0/H)*Y2(KHI))*(H*H/6.0)
  

  RETURN

END SUBROUTINE SPLINT

  
