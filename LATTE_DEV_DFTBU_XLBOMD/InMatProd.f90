subroutine InMatProd(Sum,A,B,HDIM)

IMPLICIT NONE
INTEGER, PARAMETER        :: PREC = 8
INTEGER, intent(in)  :: HDIM
REAL(PREC), intent(in)  :: A(HDIM,HDIM), B(HDIM,HDIM)
REAL(PREC), intent(out) :: Sum
INTEGER           :: I,J

Sum = 0.D0
DO J = 1,HDIM
DO I = 1,HDIM
  Sum = Sum + A(I,J)*B(I,J)
ENDDO
ENDDO

end subroutine InMatProd
