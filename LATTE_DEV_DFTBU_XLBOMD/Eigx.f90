subroutine Eigx(A,Q,ee,TYPE,HDIM)

implicit none
integer,      parameter       :: PREC = 8
integer,      intent(in)      :: HDIM
integer    :: INFO
integer    :: DIAG_LWORK
real(PREC) :: DIAG_WORK(1+6*HDIM+2*HDIM*HDIM)
integer    :: DIAG_LIWORK
integer    :: DIAG_IWORK(3+5*HDIM)
real(PREC), intent(in)      :: A(HDIM,HDIM)
real(PREC), intent(out)     :: Q(HDIM,HDIM), ee(HDIM)
real(PREC) :: EVECS(HDIM,HDIM), EVALS(HDIM)
character(LEN=1), parameter :: JOBZ = "V",  UPLO = "U"
character(1), intent(in)    :: TYPE  ! 'N'(for eigenvaules only) or 'V' (otherwise)

DIAG_LWORK = 1+6*HDIM+2*HDIM*HDIM
DIAG_LIWORK = 3 + 5*HDIM

EVECS = A
     CALL DSYEVD(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK, INFO)
Q = EVECS
ee = EVALS

end subroutine Eigx

