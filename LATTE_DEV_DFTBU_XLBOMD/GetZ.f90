subroutine GetZ(Z,S,HDIM)

implicit none
integer,    parameter       :: PREC = 8
integer,    intent(in)      :: HDIM
integer                     :: IGENX
real(PREC), intent(in)      :: S(HDIM,HDIM)
real(PREC), intent(out)     :: Z(HDIM,HDIM)
real(PREC)                  :: UMAT(HDIM,HDIM)

real(PREC)    :: NONO_EVALS(HDIM), NONO_WORK(1 + 6*HDIM + 2*HDIM*HDIM) 
real(PREC)    :: NONOTMP(HDIM,HDIM)
integer       :: NONO_IWORK(3+5*HDIM)
integer       :: I, J, INFO, NONO_LWORK
real(PREC)    :: INVSQRT,ERR_CHECK,NUMTHRESH

NONO_LWORK = 1 + 6*HDIM + 2*HDIM*HDIM

UMAT = S
if (PREC.eq.4) then
  call SSYEV("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
         NONO_LWORK, INFO)
else
  call DSYEV("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
         NONO_LWORK, INFO)
endif

do I = 1, HDIM
  INVSQRT = 1/sqrt(NONO_EVALS(I))
  do J = 1, HDIM
     NONOTMP(J,I) = UMAT(J,I)*INVSQRT
  enddo
enddo

if (PREC.eq.4) then
  call SGEMM('N', 'T', HDIM, HDIM, HDIM, 1.0, &
         NONOTMP, HDIM, UMAT, HDIM, 0.0, Z, HDIM)
else
  call DGEMM('N', 'T', HDIM, HDIM, HDIM, 1.D0, &
         NONOTMP, HDIM, UMAT, HDIM, 0.D0, Z, HDIM)
endif

end subroutine GetZ

