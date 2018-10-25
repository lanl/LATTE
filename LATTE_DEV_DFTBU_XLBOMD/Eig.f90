subroutine Eig(A,Q,ee,TYPE,HDIM)

implicit none
integer,      parameter       :: PREC = 8
integer,      intent(in)      :: HDIM
real(PREC),   intent(in)      :: A(HDIM,HDIM)
real(PREC),   intent(out)     :: Q(HDIM,HDIM), ee(HDIM)
character(1), intent(in)      :: TYPE  ! 'N'(for eigenvaules only) or 'V' (otherwise)

!real(PREC)    :: NONO_EVALS(HDIM), NONO_WORK(1 + 6*HDIM + 2*HDIM*HDIM)
real(PREC)    :: NONO_EVALS(HDIM), NONO_WORK(8 + 8*HDIM + 2*HDIM*HDIM)
!real(PREC)    :: LI_WORK(3 + 5*HDIM)
real(PREC)    :: LI_WORK(6 + 8*HDIM)
integer       :: I, J, INFO, NONO_LWORK, LIWORK

!NONO_LWORK = 1 + 6*HDIM + 2*HDIM*HDIM
!LIWORK = 1 + 6*HDIM + 2*HDIM*HDIM
NONO_LWORK = 8 + 8*HDIM + 2*HDIM*HDIM
LIWORK = 8 + 8*HDIM + 2*HDIM*HDIM

Q = A
if (PREC.eq.4) then
  call SSYEV('V', 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
         NONO_LWORK, INFO)
!  call SSYEVD('V', 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
!         NONO_LWORK, LI_WORK,LIWORK,INFO)
else
  call DSYEV('V', 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
         NONO_LWORK, INFO)
!  call DSYEVD('V', 'U', HDIM, Q, HDIM, ee, NONO_WORK, &
!         NONO_LWORK, LI_WORK,LIWORK,INFO)
endif

end subroutine Eig

