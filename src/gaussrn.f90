FUNCTION GAUSSRN(MEAN,STDDEV)

  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  !
  ! Based on GASDEV from Numerical Recipes
  !
  ! We generate 2 random numbers at a time so we have to be
  ! a little careful. 
  !

  REAL(LATTEPREC) :: MEAN, STDDEV
  REAL(LATTEPREC) :: GAUSSRN
  REAL(LATTEPREC) :: V1, V2, RN(2), R, FAC, Y1, Y2
  REAL(LATTEPREC), SAVE :: G2

  IF (SETTH .EQ. 0) THEN

     R = TWO

     DO WHILE (R .GE. ONE) 

        CALL RANDOM_NUMBER(RN)

        V1 = TWO*RN(1) - ONE
        V2 = TWO*RN(2) - ONE

        R = V1*V1 + V2*V2

     ENDDO

     FAC = SQRT( -TWO * LOG(R) / R )
     Y1 = V1 * FAC
     Y2 = V2 * FAC

     G2 = (Y1 * STDDEV) + MEAN
     GAUSSRN = (Y2 * STDDEV) + MEAN

     SETTH = 1

  ELSE 

     GAUSSRN = G2

     SETTH = 0

  ENDIF

  RETURN

END FUNCTION GAUSSRN
