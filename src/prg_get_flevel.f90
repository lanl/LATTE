!> Bisection method to get the chemical potential
!!
!!
SUBROUTINE prg_get_flevel(eigenvalues,norb,inkbt,inbndfill,tol,Ef,verbose)

  !USE CONSTANTS_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER       :: DP = 8
  INTEGER                  ::  i, j, k, m
  INTEGER, intent(in)      ::  verbose
  INTEGER                  ::  norb
  REAL(DP)                 ::  Ft1, Ft2, Kb, Prod
  REAL(DP)                 ::  T, step, tol, nel, Fermi
  REAL(DP), INTENT(in)     ::  inbndfill, eigenvalues(norb), inkbt
  REAL(DP), INTENT(inout)  ::  Ef

  IF(VERBOSE >= 1)WRITE(*,*)inkbt,inbndfill,tol,Ef,norb
  nel = inbndfill*2.0_dp*norb
  Ef=eigenvalues(1)
  step=ABS((eigenvalues(norb)-eigenvalues(1)))
  Ft1=0.0_dp
  Ft2=0.0_dp
  prod=0.0_dp

  IF(VERBOSE >= 3)WRITE(*,*)"Chemical potential =",ef
  IF(VERBOSE >= 3)WRITE(*,*)"Number of electrons =",nel

  !Sum of the occupations
  DO i=1,norb
    fermi = 1.0_dp/(1.0_dp+exp((eigenvalues(i)-ef)/(inkbt)))
     ft1 = ft1 + 2.0_dp*fermi
  ENDDO
  ft1=ft1-nel

  DO m=1,1000001

     IF(m.GT.1000000)THEN
        STOP "Bisection method in prg_get_flevel not converging ..."
     ENDIF

     IF(ABS(ft1).LT.tol)THEN !tolerance control
        RETURN
     ENDIF

     ef = ef + step

     IF(VERBOSE >= 3)WRITE(*,*)"Chemical potential =",ef

     ft2=0.0_dp

     !New sum of the occupations
     DO i=1,norb
        fermi = 1.0_dp/(1.0_dp+exp((eigenvalues(i)-ef)/(inkbt)))
        ft2 = ft2 + 2.0_dp*fermi
     ENDDO

     ft2=ft2-nel

     IF(VERBOSE >= 3)WRITE(*,*)"Occupation error =",ft2

     !Product to see the change in sign.
     prod = ft2*ft1

     IF(prod.LT.0)THEN
        ef=ef-step
        step=step/2.0_dp !If the root is inside we shorten the step.
     ELSE
        ft1=ft2  !If not, Ef moves forward.
     ENDIF

  ENDDO

END SUBROUTINE prg_get_flevel
