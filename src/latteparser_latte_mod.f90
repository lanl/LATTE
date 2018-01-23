!> LATTE parser.
!! \ingroup LATTE
!! This module is used to parse all the necessary input variables for a LATTE TB run (SCF/OPT/MD)
!! Adding a new input keyword to the parser:
!! - If the variable is real, we have to increase nkey_re.
!! - Add the keyword (character type) in the keyvector_re vector.
!! - Add a default value (real type) in the valvector_re.
!! - Define a new variable int the latte type and pass the value through valvector_re(num)
!!   where num is the position of the new keyword in the vector.
!! - Use DUMMY= as a placeholder. This variable will be ignored by not searched by the parser.
!!
module latteparser_latte_mod

  use CONSTANTS_MOD
  use SETUPARRAY
  use PPOTARRAY
  use NEBLISTARRAY
  use COULOMBARRAY
  use SPARSEARRAY
  use RELAXCOMMON
  use MDARRAY
  USE KSPACEARRAY

  use openfiles_mod
  use kernelparser_mod

#ifdef PROGRESSON
  use bml
#endif

  implicit none

  private

  integer, parameter :: dp = latteprec

  public :: parse_control, parse_md, parse_kmesh

#ifdef PROGRESSON
!> General latte input variables type.
!!
type, public :: latte_type

  !> Name of the current job.
  character(20) :: jobname

  !> Verbosity level.
  integer :: verbose

  !> Threshold values for matrix elements.
  real(dp) :: threshold

  !> Max nonzero elements per row for every row see \cite Mniszewski2015 .
  integer :: mdim

  !> Matrix format (Dense or Ellpack).
  character(20) :: bml_type

  !> Distribution mode (sequential, distributed, or graph_distributed).
  character(20) :: bml_dmode

  !> Coulomb Accuracy.
  real(dp) :: coul_acc

  !> Pulay mixing coefficient.
  real(dp) :: pulaycoeff

  !> Linear mixing coefficient.
  real(dp) :: mixcoeff

  !> Coulomb Accuracy.
  integer :: mpulay

  !> Maximum SCF iterations.
  integer :: maxscf

  !> SCF tolerance.
  real(dp) :: scftol

  !> Z Matrix calculation type.
  character(20) :: ZMat

  !> Solver method
  character(20) :: method

  !> Estimated ration between real & k space time efficiency.
  real(dp) :: timeratio

  !> Total number of steps for MD simulation.
  integer :: mdsteps

  !> Total number of steps for MD simulation.
  real(dp) :: timestep

  !> Total number of steps for MD simulation.
  character(100) :: parampath

  !> File containing coordinates.
  character(100) :: coordsfile

  !> File containing coordinates.
  integer :: nlisteach

  !> Restart calculation.
  logical :: restart

  !> Restart calculation.
  real(dp) :: efermi


end type latte_type

type(latte_type), public :: lt

#endif

contains

  !> The parser for Latte General input variables.
  !!
  subroutine parse_control(filename)

    use FERMICOMMON

    implicit none
    integer, parameter :: nkey_char = 6, nkey_int = 51, nkey_re = 21, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'JobName=','BASISTYPE=','SP2CONV=','RELAXTYPE=','PARAMPATH=','COORDSFILE=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob','NONORTHO','REL','SD','./TBparam','./bl/inputblock.dat']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'xControl=', 'DEBUGON=', 'FERMIM=', 'CGORLIB=', 'NORECS=', 'ENTROPYKIND=',&
         'PPOTON=', 'VDWON=', 'SPINON=', 'ELECTRO=', 'ELECMETH=', 'MAXSCF=',& !12
         'MINSP2ITER=','FULLQCONV=','QITER=','ORDERNMOL=','SPARSEON=','THRESHOLDON=',& !18
         'FILLINSTOP=','BLKSZ=','MSPARSE=','LCNON=','LCNITER=','RELAX=','MAXITER=',& !25
         'MDON=','PBCON=','RESTART=','CHARGE=','XBO=','XBODISON=','XBODISORDER=','NGPU=',& !33
         'KON=','COMPFORCE=','DOSFIT=','INTS2FIT=','NFITSTEP=','QFIT=',& !39
         'PPFITON=','ALLFITON=','PPSTEP=','BISTEP=','PP2FIT=','BINT2FIT=','PPNMOL=',& !46
         'PPNGEOM=','PARREP=','VERBOSE=','MIXER=','RESTARTLIB=']
    integer :: valvector_int(nkey_int) = (/ &
         1,0,6,1,1,1, &
         1,0,0,1,0,250, &
         22,0,1,0,0,1, &
         100,4,3000,0,4,0,100, &
         1,1,0,0,1,1,5,2, &
         0,1,0,1,5000,0,&
         0,0,500,500,2,6,10,&
         200,0,0,0,0 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'CGTOL=','KBT=','SPINTOL=','ELEC_ETOL=','ELEC_QTOL=','COULACC=','COULCUT=', 'COULR1=',& !8
         'BREAKTOL=','QMIX=','SPINMIX=','MDMIX=','NUMTHRESH=','CHTOL=','SKIN=',& !15
         'RLXFTOL=','BETA=','MCSIGMA=','PPBETA=','PPSIGMA=','ER='] !21
    real(dp) :: valvector_re(nkey_re) = (/&
         1.0e-6,0.0,1.0e-4,0.001,1.0e-8,1.0e-6,-500.0, 500.0,&
         1.0e-6,0.25,0.25,0.25,1.0e-6,0.01,1.0,&
         1.0e-7,1000.0,0.2,1000.0,0.01,1.0/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'LIBINIT=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'CONTROL{', '}']

    call parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)


    JOB = valvector_char(1)

    ! CONTROL determines how the density matrix is going to be
    ! calculated: 1 = diagonalization, 2 = SP2 purification,
    ! 3 = recursive expansion of the Fermi operator, 4 = SP2T,
    ! 5 = SP2 Fermi
    !

    CONTROL = valvector_int(1)

    !
    ! BASISTYPE can equal "ORTHO" OR "NONORTHO",
    !

    BASISTYPE = valvector_char(2)

    if (BASISTYPE .ne. "ORTHO" .and. BASISTYPE .ne. "NONORTHO") then
       print*, "Error defining basis type (ortho/nonortho)"
       stop
    endif

    DEBUGON = valvector_int(2)


    !
    ! Read the order of the recursion in the expansion of the Fermi
    ! operator, M.
    !

    FERMIM = valvector_int(3)

    ! If we're using the expansion of the Fermi operator, we can
    ! use a LAPACK routine or Niklasson's conjugate gradient method to
    ! solve AX = B. CGORLIB: 0 = LAPACK, 1 = conjugate gradient
    ! CGTOL = the user-supplied tolerance for the CG solution of AX = B

    CGORLIB = valvector_int(4); CGTOL = valvector_re(1)

    CGTOL2 = CGTOL*CGTOL

    ! Electronic temperature, in eV

    KBT = valvector_re(2)

    !
    ! Read the number of recursions for the truncated, finite
    ! temperature SP2 algorithm
    !

    NORECS = valvector_int(5)

    !
    ! What kind of entropy are we going to use in a finite Te calculation
    !
    ! ENTROPYKIND = 0 : none
    ! ENTROPYKIND = 1 : exact for Fermi-Dirac occupation
    ! ENTROPYKIND = 2 : Different form of exact expression that may be useful
    ! when using CONTROL = 5
    ! ENTROPYKIND = 3 : 4th order expansion of exact form (no diag)
    ! ENTROPYKIND = 4 : 8th order expansion of exact form (no diag)
    !

    ENTROPYKIND = valvector_int(6)

    !
    ! Do we want long-range C/R^6 tails?
    !
    ! PPOTON = 1: Turn on pairwise interaction
    ! PPOTON = 0: Turn it off (useful for fitting)
    !
    ! VDWON = 0: No C/R^6 tails
    ! VDWON = 1: Use tails
    !

    PPOTON = valvector_int(7);  VDWON = valvector_int(8)


    !
    ! Are we doing a spin-polarized calculation?
    ! SPINON = 1 = yes
    ! SPINON = 0 = no

    SPINON = valvector_int(9); SPINTOL = valvector_re(3)

    !
    ! Controls for electrostatics:
    !
    ! ELECTRO: 0 = LCN is applied, 1 = charge dependent TB on
    ! ELECMETH: 0 = Ewald  summation, 1 = All real space
    ! ELEC_ETOL: Tolerance on energy when determining charges (not implemented)
    ! ELEC_QTOL: Tolerance on charges during self-consistent calc
    !

    ELECTRO = valvector_int(10); ELECMETH = valvector_int(11)
    ELEC_ETOL = valvector_re(4); ELEC_QTOL = valvector_re(5)

    !
    ! COULACC: Accuracy for the Ewald method (1.0e-4 works)
    ! COULCUT: If we're using the Ewald method, this is the cut-off for the
    ! real space part. If we're doing it all in real space, this is the radial
    ! cut-off for the sum.
    ! COULR1: If we're doing it in real space, the cut-off tail on 1/R is
    ! applied here at terminated at COULCUT.
    !

    COULACC = valvector_re(6); COULCUT = valvector_re(7); COULR1 = valvector_re(8)

    !
    ! MAXSCF:  Maximum number of SCF cycles
    !

    MAXSCF = valvector_int(12)

    !
    ! BREAKTOL: Tolerance for breaking SP2 loops
    ! MINSP2ITER: Minimum number of iterations during SP2 purification
    !

    BREAKTOL = valvector_re(9); MINSP2ITER = valvector_int(13); SP2CONV = valvector_char(3)

    !
    ! FULLQCONV: 0 = We'll run QITER SCF cycles during MD, = 1, we'll run
    ! SCF cycles until we've reached ELEC_QTOL. Only important for MD
    ! QITER: Number of SCF cycles we're going to run at each MD time step
    !

    FULLQCONV = valvector_int(14); QITER = valvector_int(15)

    !
    ! QMIX AND SPINMIX are the coefficients for the linear mixing of
    ! new and old charge and spin densities, respectively, during SCF cycles
    !

    QMIX = valvector_re(10); SPINMIX = valvector_re(11); MDMIX = valvector_re(12)

    !
    ! ORDERNMOL: Turn on molecule-ID-based density matrix blocking
    !

    ORDERNMOL = valvector_int(16)

    !
    ! SPARSEON: 0 = all dense matrix stuff, 1 = use CSR format and
    ! Gustavson's algorithm for matrix-matrix multiplication
    ! THRESHOLDON: 0 = do not throw away elements; 1 = throw away elements
    ! NUMTHRESH: If THRESHOLDON = 1 throw away element whose absolute value is
    ! smaller than NUMTHRESH
    ! FILLINSTOP: Number of purification cycles beyond which we stop allowing
    ! for further fill-in
    !

    SPARSEON = valvector_int(17); THRESHOLDON = valvector_int(18); NUMTHRESH = valvector_re(13)

#ifdef PROGRESSON

    if(SPARSEON == 0)then
      lt%bml_type = bml_matrix_dense
    elseif(SPARSEON == 1)then
      lt%bml_type = bml_matrix_ellpack
    else
      STOP 'SPARSEON > 1 yet not implemented'
    endif

    if(THRESHOLDON == 0)then
      lt%threshold = 0.0_dp
    elseif(THRESHOLDON == 1)then
      lt%threshold = NUMTHRESH
    else
      STOP 'THRESHOLDON > 1 yet not implemented'
    endif

#endif

    FILLINSTOP = valvector_int(19); BLKSZ = valvector_int(20)

    !
    ! MSPARSE: value for M when SPARSEON = 1, used by sp2 sparse algorithm
    !          0 = value for M is not known, defaults to N
    !

    MSPARSE = valvector_int(21)

#ifdef PROGRESSON

    if(MSPARSE == 0)then
      lt%mdim = -1  !Defaults to N
    elseif(MSPARSE > 0)then
      lt%mdim = MSPARSE
    else
      STOP 'MSPARSE cannot be negative'
    endif

#endif

    !
    ! LCNON: 0 = during charge neutral MD simulations we'll run LCNITER SCF
    ! cycles at each time step, 1 = we'll run SCF cycles until CHTOL is reached
    ! LCNITER: Number of SCF cycles to achieve LCN at each MD time step
    ! CHTOL: Tolerance on atomic charges (Mulliken) before LCN is declared
    !

    LCNON = valvector_int(22); LCNITER = valvector_int(23); CHTOL = valvector_re(14)

    !
    ! Read the SKIN for the neighbor list (Angstrom)
    !

    SKIN = valvector_re(15)

    !
    ! RELAXME: 0 = Don't run relaxation, 1 = relax geometry
    ! RELTYPE: SD = steepest descent, CG = conjugate gradient
    ! MXRLX: Maximum number of steps in the geometry optimization
    ! RLXFTOT: Run optimization until all forces are less than RLXFTOL
    !

    RELAXME = valvector_int(24); RELTYPE = valvector_char(4) ;MXRLX = valvector_int(25)
    RLXFTOL = valvector_re(16)

    !
    ! MDON: 0 = Molecular dynamics off, 1 = Molecular dynamics on
    ! (MD is controlled using the file MDcontroller)
    !

    MDON = valvector_int(26)

    !
    ! PBCON: 1 = full periodic boundary conditions, 0 = gas phase: no pbc and
    ! electrostatics done all in real space
    !

    PBCON = valvector_int(27)

    RESTART = valvector_int(28)

    ! Add or remove electrons. 2+ -> charge = +2 since TOTNE = TOTNE - CHARGE

    CHARGE = valvector_int(29)

    !
    ! XBOON: 0 = Niklasson's extended Lagrangian Born-Oppenheimer MD off,
    ! 1 = on.
    !

    XBOON = valvector_int(30)

    !
    ! XBODISON: We have the option of turning on damping for the XBO
    ! to remedy the accumulation of noise. 0 = off, 1 = on.
    !

    XBODISON = valvector_int(31)

    !
    ! XBODISORDER: = Order of the damping function (1 - 9)
    !

    XBODISORDER = valvector_int(32)

    FITON = 0

    !
    ! Read in the number of GPUs per node
    !

    NGPU = valvector_int(33)

    ! Are we doing k-space?

    KON = valvector_int(34)

    ! Do we want to calculate forces too (not always necessary when fitting)

    COMPFORCE = valvector_int(35)

    ! Turn on the simulated annealing subroutine to fit DOS


    DOSFITON = valvector_int(36); INT2FIT = valvector_int(37)
    MCBETA = valvector_re(17); NFITSTEP = valvector_int(38); QFIT = valvector_int(39)
    MCSIGMA = valvector_re(18)

    PPFITON = valvector_int(40)

    ALLFITON = valvector_int(41)

    PPNFITSTEP = valvector_int(42); BINFITSTEP = valvector_int(43)
    PP2FIT = valvector_int(44); BINT2FIT = valvector_int(45)

    PPBETA = valvector_re(19); PPSIGMA = valvector_re(20)
    PPNMOL = valvector_int(46); PPNGEOM = valvector_int(47)

    PARREP = valvector_int(48)

    ! Verbosity level to control general output

    VERBOSE = valvector_int(49)

    ! If Pulay Mixer

    MIXER = valvector_int(50)

    ! Restart option for latte lib.

    MIXER = valvector_int(51)

    ! Dielectric constant

    RELPERM = valvector_re(21)

    ! Coordinates and parameter paths

    PARAMPATH = valvector_char(5)
    COORDSFILE = valvector_char(6)

  end subroutine parse_control


  !> The parser for Latte General input variables.
  !!
  subroutine parse_md(filename)

    implicit none
    integer, parameter :: nkey_char = 3, nkey_int = 18, nkey_re = 10, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'RNDIST=','SEEDINIT=','NPTTYPE=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'GAUSSIAN','UNIFORM','ISO']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'MAXITER=', 'UDNEIGH=', 'DUMPFREQ=','RSFREQ=', 'WRTFREQ=', 'TOINITTEMP5=', 'THERMPER=',& !7
         'THERMRUN=', 'NVTON=', 'NPTON=', 'AVEPER=', 'SEED=', 'SHOCKON=',&
         'SHOCKSTART=','SHOCKDIR=','MDADAPT=','GETHUG=','RSLEVEL=']
    integer :: valvector_int(nkey_int) = (/ &
         5000,1,250,500,25,1,500, &
         50000,0,0,1000,54,0, &
         100000,1,0,0,0/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'DT=','TEMPERATURE=','FRICTION=','PTARGET=','UPARTICLE=','USHOCK=','C0=', 'E0=',&
         'V0=','P0=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.25,300.00,1000.0,0.0,500.0,-4590.0,1300.0,-795.725,&
         896.984864,0.083149/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'DUMMY=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'MDCONTROL{', '}']

    call parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)


    !
    ! MAXITER = run this many MD time steps
    !

    MAXITER = valvector_int(1)

    !
    ! UDNEIGH = update the neighbor lists every UDNEIGH time steps
    !

    UDNEIGH = valvector_int(2)
    !   write(*,*)"UDNEIGH",UDNEIGH
    !
    ! DT = size of the time step in fs
    !

    DT = valvector_re(1)
    !   write(*,*)"DT",DT
    !
    ! TTARGET = temperature in K were initialize and aim for during NVT MD
    ! RNDIST = Type of distribution of random numbers used to initialize T
    !        = GAUSSIAN or UNIFORM
    ! SEEDINIT = Type of seed used in the generation of random numbers
    !          = RANDOM - seed changes every time
    !          = DEFAULT - use the same, default seed every time
    !

    TTARGET = valvector_re(2); RNDIST = valvector_char(1); SEEDINIT = valvector_char(2)
    !   write(*,*)"TTARGET,RNDIST,SEEDINIT",TTARGET,RNDIST,SEEDINIT
    !
    ! DUMPFREQ: Write a dump file every DUMPFREQ time steps
    !

    DUMPFREQ = valvector_int(3)
    !   write(*,*)"DUMPFREQ",DUMPFREQ
    !
    ! RSFREQ: Write a restart file every RSFREQ time steps
    !

    RSFREQ = valvector_int(4)
    !   write(*,*)"RSFREQ",RSFREQ
    !
    ! WRTFREQ: Output energy and temperature every WRTFREQ time steps
    !

    WRTFREQ = valvector_int(5)
    !   write(*,*)"WRTFREQ",WRTFREQ
    !
    ! TOINITTEMP: Whether or not we are going to initialize velocities
    ! using a random number generator (sometimes during a restart we
    ! may not want to reinitialize the temperature
    !

    TOINITTEMP = valvector_int(6)
    !   write(*,*)"TOINITTEMP",TOINITTEMP
    !
    ! THERMPER: If we're running NVT, rescale velocities every THERMPER
    ! time steps.
    !

    THERMPER = valvector_int(7)
    !   write(*,*)"THERMPER",THERMPER
    !
    ! THERMRUN: Thermalize over this many time steps when NVT is on
    !

    THERMRUN = valvector_int(8)
    !   write(*,*)"THERMRUN",THERMRUN
    !
    ! NVTON: 0 = running NVE MD, 1 = running NVT MD
    ! AVEPER: Average the temperature over AVEPER time steps when determining
    ! how to rescale velocities
    !

    NVTON = valvector_int(9); NPTON = valvector_int(10)
    !   write(*,*)"NVTON,NPTON",NVTON,NPTON
    AVEPER = valvector_int(11); FRICTION = valvector_re(3); SEEDTH = valvector_int(12)
    !   write(*,*)"AVEPER,FRICTION,SEEDTH",AVEPER,FRICTION,SEEDTH


    if (NVTON .eq. 1 .and. NPTON .eq. 1) then
       write(6,*) "You can't have NVTON = 1 and NPTON = 1"
       write(6,*) "STOP!!!"
       stop
    endif

    ! PTARGET = Target pressure (in GPa) when running NPT
    ! NPTTYPE = ISO or ANISO

    PTARGET = valvector_re(4); NPTTYPE = valvector_char(3)
    !   write(*,*)"PTARGET,NPTTYPE",PTARGET,NPTTYPE

    !
    ! The following are for the Hugoniostat
    !

    ! On (1) or off (0)?

    SHOCKON = valvector_int(13)
    !   write(*,*)"SHOCKON",SHOCKON
    !
    ! SHOCKSTART = the MD iteration where we will start to compress
    ! the iteration when we stop depends on the size of the block and Us
    !

    SHOCKSTART = valvector_int(14)
    !   write(*,*)"SHOCKSTART",SHOCKSTART
    !
    ! SHOCKDIR is the cartensian direction (1 = X, 2 = Y, 3 = Z),
    ! parallel to which we're going to compress uniaxially
    !

    SHOCKDIR = valvector_int(15)
    !   write(*,*)"SHOCKDIR",SHOCKDIR
    !
    ! And finally, the particle and shock velocities
    ! IN UNITS OF METRES PER SECOND
    !

    UPARTICLE = valvector_re(5); USHOCK = valvector_re(6); C0 = valvector_re(7)
    !   write(*,*)"UPARTICLE,USHOCK,C0",UPARTICLE,USHOCK,C0
    ! Adapt SCF on the fly?

    MDADAPT = valvector_int(16)
    !   write(*,*)"MDADAPT",MDADAPT
    ! Calculating Hugoniot points?

    GETHUG = valvector_int(17)

    RSLEVEL = valvector_int(18)

    !   write(*,*)"GETHUG",GETHUG
    E0 = valvector_re(8); V0 = valvector_re(9); P0 = valvector_re(10)
    !   write(*,*)"E0,V0,P0",E0,V0,P0

  end subroutine parse_md


  !> The parser for K Mesh input variables.
  !!
  subroutine parse_kmesh(filename)

    implicit none
    integer, parameter :: nkey_char = 1, nkey_int = 3, nkey_re = 3, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'DUMMY=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'DUMMY']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'NKX=','NKY=','NKZ=']
    integer :: valvector_int(nkey_int) = (/ &
         1,1,1/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'KSHIFTX=','KSHIFTY=','KSHIFTZ=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0,0.0,0.0/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'DUMMY=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'KMESH{', '}']

    call parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

      NKX= valvector_int(1); NKY= valvector_int(2); NKZ=valvector_int(3)
      KSHIFT(1)= valvector_re(1); KSHIFT(2)= valvector_re(2); KSHIFT(3)= valvector_re(3)


  end subroutine parse_kmesh


end module latteparser_latte_mod
