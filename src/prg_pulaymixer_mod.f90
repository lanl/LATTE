!> Pulay mixer mode. From the progress library
!! \ingroup PROGRESS
!! Gets the best coefficient for mixing the charges during scf.
!! \todo add the density matrix mixer.
module prg_pulaymixer_mod

  use kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  type, public :: mx_type

    !> Type or mixing scheme to be used (Linear or Pulay)
    character(20) :: mixertype

    !> Verbosity level
    integer       :: verbose

    !> Pulay dimension for matrix
    integer       :: mpulay

    !> Coefficient for mixing
    real(dp)      :: mixcoeff

    !> Mixer on or off (Not implemented)
    logical       :: mixeron

  end type mx_type

  public :: prg_parse_mixer, prg_qmixer, prg_linearmixer

contains

  !> The parser for the mixer routines.
  !!
  subroutine prg_parse_mixer(input,filename)
    implicit none
    type(mx_type), intent(inout) :: input
    integer, parameter :: nkey_char = 1, nkey_int = 2, nkey_re = 1, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
      'MixerType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
      'Linear']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
    'Verbose=','MPulay=']
    integer :: valvector_int(nkey_int) = (/ &
       0,5/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
      'MixCoeff=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.25 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
      'MixerON=']
    logical :: valvector_log(nkey_log) = (/&
     .true./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
      'MIXER{', '}']

    call parsing_kernel(keyvector_char,valvector_char&
    ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
    keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    input%mixertype = valvector_char(1)

    !Integers
    input%verbose = valvector_int(1)
    input%mpulay = valvector_int(2)

    !Logicals
    input%mixeron = valvector_log(1)

    !Reals
    input%mixcoeff = valvector_re(1)

  end subroutine prg_parse_mixer


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
  subroutine prg_qmixer(charges,oldcharges,dqin,dqout,scferror,piter,pulaycoef,mpulay,verbose)

    implicit none
    integer :: i,j,info,s,k,n
    real(dp) :: alpha,coeff
    real(dp), intent(in) :: pulaycoef
    real(dp), intent(inout) :: scferror
    real(dp) :: error, errop
    integer :: piter
    real(dp) :: dt, atenuacion
    integer :: nint, na, nb
    integer, intent(in) :: mpulay,verbose
    real(dp), allocatable :: d(:),dl(:),dnewin(:)
    real(dp), allocatable :: dnewout(:)
    real(dp), intent(inout) :: charges(:)
    real(dp), allocatable, intent(inout) :: oldcharges(:)
    real(dp), allocatable, intent(inout) :: dqin(:,:),dqout(:,:)
    real(dp), allocatable :: coef(:,:),b(:),ipiv(:)

    n=size(charges)

    alpha = pulaycoef !the coefficient for mixing

    if(allocated(oldcharges).eqv..false.)then
      allocate(oldcharges(n),dqin(n,mpulay),dqout(n,mpulay))
    endif

    if(allocated(dqin).eqv..false.)then
      allocate(dqin(n,mpulay),dqout(n,mpulay))
    endif

    s=min(piter-1,mpulay) !mpulay is the iteration number

    if(piter.eq.1) then
      charges=(1.0_dp-alpha)*oldcharges + alpha*charges
      scferror = maxval(abs(charges-oldcharges))
      if(verbose.ge.1)then
        write(*,*)"SCF error =", scferror
      endif
       oldcharges=charges
    else

      allocate(d(n),dnewin(n),dnewout(n))

      d=charges

      allocate(coef(s+1,s+1)) !Allocating the coeffs matrix
      allocate(b(s+1))
      allocate(ipiv(s+1))

      if(piter.le.mpulay+1)then  !If piter=6 => mpulay=5
        dqin(:,piter-1)=oldcharges(:)
        dqout(:,piter-1)=d(:)
      endif

      if(piter.gt.mpulay+1)then

        do j=1,s-1
          dqin(:,j)=dqin(:,j+1)
          dqout(:,j)=dqout(:,j+1)
        enddo

        dqin(:,s)=oldcharges(:)
        dqout(:,s)=d(:)

      endif

      coef=0.0_dp

      do i=1,s+1
        coef(s+1,i)=-1.0d0
        coef(i,s+1)=-1.0d0
        b(i)=0
      enddo
      b(s+1)=-1.0d0
      coef(s+1,s+1)=0.0_dp

      do i=1,s
        do j=1,s
          do k=1,n
            coef(i,j)=coef(i,j)+(dqout(k,i)-dqin(k,i))*(dqout(k,j)-dqin(k,j))
          enddo
        enddo
      enddo

      if(verbose.ge.1)then
        write(*,*)"coefs"
        do i=1,s+1
          write(*,'(10f12.5)')(coef(i,j),j=1,s+1)
        enddo
        write(*,*)"dqin"
        write(*,'(10f12.5)')(dqin(n,j),j=1,s)
      endif

      call dgesv(s+1,1,coef,s+1,ipiv,b,s+1,info)

      if(info.ne.0) stop 'singular matrix in pulay'

      dnewin=0.0_dp
      dnewout=0.0_dp

      if(verbose.ge.1)then
        write(*,*)"eigen coefs"
        write(*,'(6f10.5)')(b(j),j=1,s)
      endif

      do j=1,s
        dnewin(:)=dnewin(:)+b(j)*dqin(:,j)
        dnewout(:)= dnewout(:)+b(j)*dqout(:,j)
      enddo

      d=(1.0_dp-alpha)*dnewin + alpha*dnewout

      scferror = maxval(abs(d-oldcharges))

      if(verbose.ge.1)then
        write(*,*)"SCF error =", scferror
      endif

      charges=d

      oldcharges=d

    endif

  end subroutine prg_qmixer

  !> Routine to perform linear mixing.
  !! \param charges Actual charges of the system.
  !! \param oldcharges Previous scf charges.
  !! \param scferror SCF error.
  !! \param linmixcoef Mixing coefficient.
  !! \param verbose Verbosity level.
  subroutine prg_linearmixer(charges,oldcharges,scferror,linmixcoef,verbose)
    implicit none
    real(dp), intent(in) :: linmixcoef
    real(dp), intent(inout) :: scferror
    integer, intent(in) :: verbose
    real(dp), allocatable, intent(inout) :: charges(:),oldcharges(:)

    scferror = norm2(charges(:)-oldcharges(:))

    if(verbose.ge.1)then
      write(*,*)"SCF error =", scferror
    endif

    charges = (1.0_dp - linmixcoef)*oldcharges + linmixcoef*charges
    oldcharges = charges

  end subroutine prg_linearmixer

end module prg_pulaymixer_mod
