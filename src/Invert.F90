subroutine Invert(A,AI,N)

implicit none
integer, parameter                          :: PREC = 8
integer, intent(in)                         :: N
real(PREC),                 intent(in)      :: A(N,N)
real(PREC), intent(out)                     :: AI(N,N)
real(PREC)                                  :: WORK(N+N*N), C(N,N)
integer                                     :: LDA, LWORK, M, INFO, IPIV(N)
integer                                     :: I,J,K

external DGETRF
external DGETRI

AI = A
LDA = N
M = N
LWORK = N+N*N

call DGETRF(M, N, AI, LDA, IPIV, INFO)
call DGETRI(N, AI, N, IPIV, WORK, LWORK, INFO)

end subroutine Invert

