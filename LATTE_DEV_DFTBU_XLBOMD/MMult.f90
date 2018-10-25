subroutine MMult(alpha,A,B,beta,C,TA,TB,HDIM)

implicit none
integer,      parameter     :: PREC = 8
integer,      intent(in)    :: HDIM
real(PREC),   intent(inout)    :: A(HDIM,HDIM), B(HDIM,HDIM), alpha, beta
real(PREC),   intent(inout) :: C(HDIM,HDIM)
character(1), intent(in)    :: TA, TB
integer :: I,J

!write(*,*) ' FORE A(1,1) = ', A(1,1), A(2,2)
!DO I = 1,HDIM
!  A(I,I) = floor(A(I,I)*100000000.D0)/100000000.D0
!ENDDO
!write(*,*) ' EFTER A(1,1) = ', A(1,1), A(2,2)

if (PREC.eq.4) then
 call SGEMM(TA, TB, HDIM, HDIM, HDIM, alpha, &
             A, HDIM, B, HDIM, beta, C, HDIM)
else
  call DGEMM(TA, TB, HDIM, HDIM, HDIM, alpha, &
         A, HDIM, B, HDIM, beta, C, HDIM)

endif

end subroutine MMult

