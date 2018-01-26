!> Pulay mixer mode. From the progress library
!! \ingroup PROGRESS
!! Gets the best coefficient for mixing the charges during scf.
!! \todo add the density matrix mixer.
MODULE prg_pulaymixer_mod

  USE kernelparser_mod

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  TYPE, PUBLIC :: mx_type

     !> Type or mixing scheme to be used (Linear or Pulay)
     CHARACTER(20) :: mixertype

     !> Verbosity level
     INTEGER       :: verbose

     !> Pulay dimension for matrix
     INTEGER       :: mpulay

     !> Coefficient for mixing
     REAL(dp)      :: mixcoeff

     !> Mixer on or off (Not implemented)
     LOGICAL       :: mixeron

  END TYPE mx_type

  PUBLIC :: prg_parse_mixer, prg_qmixer, prg_linearmixer

CONTAINS

  !> The parser for the mixer routines.
  !!
  SUBROUTINE prg_parse_mixer(input,filename)
    IMPLICIT NONE
    TYPE(mx_type), INTENT(inout) :: input
    INTEGER, PARAMETER :: nkey_char = 1, nkey_int = 2, nkey_re = 1, nkey_log = 1
    CHARACTER(len=*) :: filename

    !Library of keywords with the respective defaults.
    CHARACTER(len=50), PARAMETER :: keyvector_char(nkey_char) = [CHARACTER(len=100) :: &
         'MixerType=']
    CHARACTER(len=100) :: valvector_char(nkey_char) = [CHARACTER(len=100) :: &
         'Linear']

    CHARACTER(len=50), PARAMETER :: keyvector_int(nkey_int) = [CHARACTER(len=50) :: &
         'Verbose=','MPulay=']
    INTEGER :: valvector_int(nkey_int) = (/ &
         0,5/)

    CHARACTER(len=50), PARAMETER :: keyvector_re(nkey_re) = [CHARACTER(len=50) :: &
         'MixCoeff=']
    REAL(dp) :: valvector_re(nkey_re) = (/&
         0.25 /)

    CHARACTER(len=50), PARAMETER :: keyvector_log(nkey_log) = [CHARACTER(len=100) :: &
         'MixerON=']
    LOGICAL :: valvector_log(nkey_log) = (/&
         .TRUE./)

    !Start and stop characters
    CHARACTER(len=50), PARAMETER :: startstop(2) = [CHARACTER(len=50) :: &
         'MIXER{', '}']

    CALL parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,TRIM(filename),startstop)

    !Characters
    input%mixertype = valvector_char(1)

    !Integers
    input%verbose = valvector_int(1)
    input%mpulay = valvector_int(2)

    !Logicals
    input%mixeron = valvector_log(1)

    !Reals
    input%mixcoeff = valvector_re(1)

  END SUBROUTINE prg_parse_mixer


  !> Mixing the charges to acelerate scf convergence.
  !! \param charges System charges.
  !! \param oldcharges Old charges of the system.
  !! \param dqin Matrix for charges history in.
  !! \param dqout Matrix for charges history out.
  !! \param scferror SCF error.
  !! \param piter scf iteration number.
  !! \param pulaycoef Coefficient for pulay mixing (generally between 0.01 and 0.1).
  !! \param mpulay Number of matrices stored (generally 3-5).
  !! \param verbose Different levels of verbosity.
  SUBROUTINE prg_qmixer(charges,oldcharges,dqin,dqout,scferror,piter,pulaycoef,mpulay,verbose)

    IMPLICIT NONE
    INTEGER :: i,j,info,s,k,n
    REAL(dp) :: alpha,coeff
    REAL(dp), INTENT(in) :: pulaycoef
    REAL(dp), INTENT(inout) :: scferror
    REAL(dp) :: error, errop
    INTEGER :: piter
    REAL(dp) :: dt, atenuacion
    INTEGER :: nint, na, nb
    INTEGER, INTENT(in) :: mpulay,verbose
    REAL(dp), ALLOCATABLE :: d(:),dl(:),dnewin(:)
    REAL(dp), ALLOCATABLE :: dnewout(:)
    REAL(dp), INTENT(inout) :: charges(:)
    REAL(dp), ALLOCATABLE, INTENT(inout) :: oldcharges(:)
    REAL(dp), ALLOCATABLE, INTENT(inout) :: dqin(:,:),dqout(:,:)
    REAL(dp), ALLOCATABLE :: coef(:,:),b(:),ipiv(:)

    n=SIZE(charges)

    alpha = pulaycoef !the coefficient for mixing

    IF(ALLOCATED(oldcharges).EQV..FALSE.)THEN
       ALLOCATE(oldcharges(n),dqin(n,mpulay),dqout(n,mpulay))
    ENDIF

    IF(ALLOCATED(dqin).EQV..FALSE.)THEN
       ALLOCATE(dqin(n,mpulay),dqout(n,mpulay))
    ENDIF

    s=MIN(piter-1,mpulay) !mpulay is the iteration number

    IF(piter.EQ.1) THEN
       charges=(1.0_dp-alpha)*oldcharges + alpha*charges
       scferror = MAXVAL(ABS(charges-oldcharges))
       IF(verbose.GE.1)THEN
          WRITE(*,*)"SCF error =", scferror
       ENDIF
       oldcharges=charges
    ELSE

       ALLOCATE(d(n),dnewin(n),dnewout(n))

       d=charges

       ALLOCATE(coef(s+1,s+1)) !Allocating the coeffs matrix
       ALLOCATE(b(s+1))
       ALLOCATE(ipiv(s+1))

       IF(piter.LE.mpulay+1)THEN  !If piter=6 => mpulay=5
          dqin(:,piter-1)=oldcharges(:)
          dqout(:,piter-1)=d(:)
       ENDIF

       IF(piter.GT.mpulay+1)THEN

          DO j=1,s-1
             dqin(:,j)=dqin(:,j+1)
             dqout(:,j)=dqout(:,j+1)
          ENDDO

          dqin(:,s)=oldcharges(:)
          dqout(:,s)=d(:)

       ENDIF

       coef=0.0_dp

       DO i=1,s+1
          coef(s+1,i)=-1.0d0
          coef(i,s+1)=-1.0d0
          b(i)=0
       ENDDO
       b(s+1)=-1.0d0
       coef(s+1,s+1)=0.0_dp

       DO i=1,s
          DO j=1,s
             DO k=1,n
                coef(i,j)=coef(i,j)+(dqout(k,i)-dqin(k,i))*(dqout(k,j)-dqin(k,j))
             ENDDO
          ENDDO
       ENDDO

       IF(verbose.GE.1)THEN
          WRITE(*,*)"coefs"
          DO i=1,s+1
             WRITE(*,'(10f12.5)')(coef(i,j),j=1,s+1)
          ENDDO
          WRITE(*,*)"dqin"
          WRITE(*,'(10f12.5)')(dqin(n,j),j=1,s)
       ENDIF

       CALL dgesv(s+1,1,coef,s+1,ipiv,b,s+1,info)

       IF(info.NE.0) STOP 'singular matrix in pulay'

       dnewin=0.0_dp
       dnewout=0.0_dp

       IF(verbose.GE.1)THEN
          WRITE(*,*)"eigen coefs"
          WRITE(*,'(6f10.5)')(b(j),j=1,s)
       ENDIF

       DO j=1,s
          dnewin(:)=dnewin(:)+b(j)*dqin(:,j)
          dnewout(:)= dnewout(:)+b(j)*dqout(:,j)
       ENDDO

       d=(1.0_dp-alpha)*dnewin + alpha*dnewout

       scferror = MAXVAL(ABS(d-oldcharges))

       IF(verbose.GE.1)THEN
          WRITE(*,*)"SCF error =", scferror
       ENDIF

       charges=d

       oldcharges=d

    ENDIF

  END SUBROUTINE prg_qmixer

  !> Routine to perform linear mixing.
  !! \param charges Actual charges of the system.
  !! \param oldcharges Previous scf charges.
  !! \param scferror SCF error.
  !! \param linmixcoef Mixing coefficient.
  !! \param verbose Verbosity level.
  SUBROUTINE prg_linearmixer(charges,oldcharges,scferror,linmixcoef,verbose)
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: linmixcoef
    REAL(dp), INTENT(inout) :: scferror
    INTEGER, INTENT(in) :: verbose
    REAL(dp), ALLOCATABLE, INTENT(inout) :: charges(:),oldcharges(:)

    scferror = NORM2(charges(:)-oldcharges(:))

    IF(verbose.GE.1)THEN
       WRITE(*,*)"SCF error =", scferror
    ENDIF

    charges = (1.0_dp - linmixcoef)*oldcharges + linmixcoef*charges
    oldcharges = charges

  END SUBROUTINE prg_linearmixer

END MODULE prg_pulaymixer_mod
