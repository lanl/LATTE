!> LATTE parser.
!! \ingroup LATTE
!! \brief This module is used to parse all the necessary input variables for a LATTE TB run (SCF/OPT/MD)
!! Adding a new input keyword to the parser:
!! - If the variable is real, we have to increase nkey_re.
!! - Add the keyword (character type) in the keyvector_re vector.
!! - Add a default value (real type) in the valvector_re.
!! - Define a new variable int the latte type and pass the value through valvector_re(num)
!!   where num is the position of the new keyword in the vector.
!! - Use DUMMY= as a placeholder. This variable will be ignored by not searched by the parser.
!!
MODULE LATTEPARSER_LATTE_MOD

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE COULOMBARRAY
  USE SPARSEARRAY
  USE RELAXCOMMON
  USE MDARRAY
  USE KSPACEARRAY

  USE OPENFILES_MOD
  USE KERNELPARSER_MOD

#ifdef PROGRESSON
  USE BML
#endif

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: DP = LATTEPREC

  PUBLIC :: PARSE_CONTROL, PARSE_MD, PARSE_KMESH

#ifdef PROGRESSON
  !> General latte input variables type.
  !!
  TYPE, PUBLIC :: LATTE_TYPE

     !> Name of the current job.
     CHARACTER(20) :: JOBNAME

     !> Verbosity level.
     INTEGER :: VERBOSE

     !> Threshold values for matrix elements.
     REAL(DP) :: THRESHOLD

     !> Max nonzero elements per row for every row see \cite Mniszewski2015 .
     INTEGER :: MDIM

     !> Matrix format (Dense or Ellpack).
     CHARACTER(20) :: BML_TYPE

     !> Distribution mode (sequential, distributed, or graph_distributed).
     CHARACTER(20) :: BML_DMODE

     !> Coulomb Accuracy.
     REAL(DP) :: COUL_ACC

     !> Pulay mixing coefficient.
     REAL(DP) :: PULAYCOEFF

     !> Linear mixing coefficient.
     REAL(DP) :: MIXCOEFF

     !> Coulomb Accuracy.
     INTEGER :: MPULAY

     !> Maximum SCF iterations.
     INTEGER :: MAXSCF

     !> SCF tolerance.
     REAL(DP) :: SCFTOL

     !> Z Matrix calculation type.
     CHARACTER(20) :: ZMAT

     !> Solver method
     CHARACTER(20) :: METHOD

     !> Estimated ration between real & k space time efficiency.
     REAL(DP) :: TIMERATIO

     !> Total number of steps for MD simulation.
     INTEGER :: MDSTEPS

     !> Total number of steps for MD simulation.
     REAL(DP) :: TIMESTEP

     !> Total number of steps for MD simulation.
     CHARACTER(100) :: PARAMPATH

     !> File containing coordinates.
     CHARACTER(100) :: COORDSFILE

     !> File containing coordinates.
     INTEGER :: NLISTEACH

     !> Restart calculation.
     LOGICAL :: RESTART

     !> Restart calculation.
     REAL(DP) :: EFERMI


  END TYPE LATTE_TYPE

  TYPE(LATTE_TYPE), PUBLIC :: LT

#endif

CONTAINS

  !> The parser for Latte General input variables.
  !!
  SUBROUTINE PARSE_CONTROL(FILENAME)

    USE FERMICOMMON

    IMPLICIT NONE
    INTEGER, PARAMETER :: NKEY_CHAR = 7, NKEY_INT = 53, NKEY_RE = 21, NKEY_LOG = 2
    CHARACTER(LEN=*) :: FILENAME

    !Library of keywords with the respective defaults.
    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'JOBNAME=','BASISTYPE=','SP2CONV=','RELAXTYPE=','PARAMPATH=','COORDSFILE=', &
         'SCLTYPE=' ]
    CHARACTER(LEN=100) :: VALVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'MyJob','NONORTHO','REL','SD','./TBparam','./bl/inputblock.dat', &
         'EXP']

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_INT(NKEY_INT) = [CHARACTER(LEN=50) :: &
         'XCONTROL=','DEBUGON=','FERMIM=','CGORLIB=','NORECS=','ENTROPYKIND=',&
         'PPOTON=','VDWON=','SPINON=','ELECTRO=', 'ELECMETH=','MAXSCF=',& !12
         'MINSP2ITER=','FULLQCONV=','QITER=','ORDERNMOL=','SPARSEON=','THRESHOLDON=',& !18
         'FILLINSTOP=','BLKSZ=','MSPARSE=','LCNON=','LCNITER=','RELAX=','MAXITER=',& !25
         'MDON=','PBCON=','RESTART=','CHARGE=','XBO=','XBODISON=','XBODISORDER=','NGPU=',& !33
         'KON=','COMPFORCE=','DOSFIT=','INTS2FIT=','NFITSTEP=','QFIT=',& !39
         'PPFITON=','ALLFITON=','PPSTEP=','BISTEP=','PP2FIT=','BINT2FIT=','PPNMOL=',& !46
         'PPNGEOM=','PARREP=','VERBOSE=','MIXER=','RESTARTLIB=','FREEZE=','xControl=']
    INTEGER :: VALVECTOR_INT(NKEY_INT) = (/ &
         1,0,6,1,1,1, &
         1,0,0,1,0,250, &
         22,0,1,0,0,1, &
         100,4,3000,0,4,0,100, &
         1,1,0,0,1,1,5,2, &
         0,1,0,1,5000,0,&
         0,0,500,500,2,6,10,&
         200,0,1,0,0,0,-1 /)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_RE(NKEY_RE) = [CHARACTER(LEN=50) :: &
         'CGTOL=','KBT=','SPINTOL=','ELEC_ETOL=','ELEC_QTOL=','COULACC=','COULCUT=', 'COULR1=',& !8
         'BREAKTOL=','QMIX=','SPINMIX=','MDMIX=','NUMTHRESH=','CHTOL=','SKIN=',& !15
         'RLXFTOL=','BETA=','MCSIGMA=','PPBETA=','PPSIGMA=','ER='] !21
    REAL(DP) :: VALVECTOR_RE(NKEY_RE) = (/&
         1.0e-6,0.0,1.0e-4,0.001,1.0e-8,1.0e-6,-500.0, 500.0,&
         1.0e-6,0.25,0.25,0.25,1.0e-6,0.01,1.0,&
         1.0e-7,1000.0,0.2,1000.0,0.01,1.0/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_LOG(NKEY_LOG) = [CHARACTER(LEN=100) :: &
         'LIBINIT=','STOPATMAXSCF=']
    LOGICAL :: VALVECTOR_LOG(NKEY_LOG) = (/&
         .FALSE.,.FALSE./)

    !Start and stop characters
    CHARACTER(LEN=50), PARAMETER :: STARTSTOP(2) = [CHARACTER(LEN=50) :: &
         'CONTROL{', '}']

    CALL PARSING_KERNEL(KEYVECTOR_CHAR,VALVECTOR_CHAR&
         ,KEYVECTOR_INT,VALVECTOR_INT,KEYVECTOR_RE,VALVECTOR_RE,&
         KEYVECTOR_LOG,VALVECTOR_LOG,TRIM(FILENAME),STARTSTOP)


    JOB = VALVECTOR_CHAR(1)

    ! CONTROL determines how the density matrix is going to be
    ! calculated: 1 = diagonalization, 2 = SP2 purification,
    ! 3 = recursive expansion of the Fermi operator, 4 = SP2T,
    ! 5 = SP2 Fermi
    !

    CONTROL = VALVECTOR_INT(1)
    IF (VALVECTOR_INT(53) > 0) THEN !Someone is using xControl=
        CALL ERRORS("latteparser_latte_mod","xControl= is not longer in use. Please use XCONTROL= instead.")
    ENDIF


    !
    ! BASISTYPE can equal "ORTHO" OR "NONORTHO",
    !

    BASISTYPE = VALVECTOR_CHAR(2)

    IF (BASISTYPE .NE. "ORTHO" .AND. BASISTYPE .NE. "NONORTHO") THEN
       CALL ERRORS("latteparser_latte_mod","Error defining basis type (ortho/nonortho)")
    ENDIF

    DEBUGON = VALVECTOR_INT(2)


    !
    ! Read the order of the recursion in the expansion of the Fermi
    ! operator, M.
    !

    FERMIM = VALVECTOR_INT(3)

    ! If we're using the expansion of the Fermi operator, we can
    ! use a LAPACK routine or Niklasson's conjugate gradient method to
    ! solve AX = B. CGORLIB: 0 = LAPACK, 1 = conjugate gradient
    ! CGTOL = the user-supplied tolerance for the CG solution of AX = B

    CGORLIB = VALVECTOR_INT(4); CGTOL = VALVECTOR_RE(1)

    CGTOL2 = CGTOL*CGTOL

    ! Electronic temperature, in eV

    KBT = VALVECTOR_RE(2)

    !
    ! Read the number of recursions for the truncated, finite
    ! temperature SP2 algorithm
    !

    NORECS = VALVECTOR_INT(5)

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

    ENTROPYKIND = VALVECTOR_INT(6)

    !
    ! Do we want long-range C/R^6 tails?
    !
    ! PPOTON = 1: Turn on pairwise interaction
    ! PPOTON = 0: Turn it off (useful for fitting)
    !
    ! VDWON = 0: No C/R^6 tails
    ! VDWON = 1: Use tails
    !

    PPOTON = VALVECTOR_INT(7);  VDWON = VALVECTOR_INT(8)


    !
    ! Are we doing a spin-polarized calculation?
    ! SPINON = 1 = yes
    ! SPINON = 0 = no

    SPINON = VALVECTOR_INT(9); SPINTOL = VALVECTOR_RE(3)

    !
    ! Controls for electrostatics:
    !
    ! ELECTRO: 0 = LCN is applied, 1 = charge dependent TB on
    ! ELECMETH: 0 = Ewald  summation, 1 = All real space
    ! ELEC_ETOL: Tolerance on energy when determining charges (not implemented)
    ! ELEC_QTOL: Tolerance on charges during self-consistent calc
    !

    ELECTRO = VALVECTOR_INT(10); ELECMETH = VALVECTOR_INT(11)
    ELEC_ETOL = VALVECTOR_RE(4); ELEC_QTOL = VALVECTOR_RE(5)

    !
    ! COULACC: Accuracy for the Ewald method (1.0e-4 works)
    ! COULCUT: If we're using the Ewald method, this is the cut-off for the
    ! real space part. If we're doing it all in real space, this is the radial
    ! cut-off for the sum.
    ! COULR1: If we're doing it in real space, the cut-off tail on 1/R is
    ! applied here at terminated at COULCUT.
    !

    COULACC = VALVECTOR_RE(6); COULCUT = VALVECTOR_RE(7); COULR1 = VALVECTOR_RE(8)

    !
    ! MAXSCF:  Maximum number of SCF cycles
    !

    MAXSCF = VALVECTOR_INT(12)

    !
    ! BREAKTOL: Tolerance for breaking SP2 loops
    ! MINSP2ITER: Minimum number of iterations during SP2 purification
    !

    BREAKTOL = VALVECTOR_RE(9); MINSP2ITER = VALVECTOR_INT(13); SP2CONV = VALVECTOR_CHAR(3)

    !
    ! FULLQCONV: 0 = We'll run QITER SCF cycles during MD, = 1, we'll run
    ! SCF cycles until we've reached ELEC_QTOL. Only important for MD
    ! QITER: Number of SCF cycles we're going to run at each MD time step
    !

    FULLQCONV = VALVECTOR_INT(14); QITER = VALVECTOR_INT(15)

    !
    ! QMIX AND SPINMIX are the coefficients for the linear mixing of
    ! new and old charge and spin densities, respectively, during SCF cycles
    !

    QMIX = VALVECTOR_RE(10); SPINMIX = VALVECTOR_RE(11); MDMIX = VALVECTOR_RE(12)

    !
    ! ORDERNMOL: Turn on molecule-ID-based density matrix blocking
    !

    ORDERNMOL = VALVECTOR_INT(16)

    !
    ! SPARSEON: 0 = all dense matrix stuff, 1 = use CSR format and
    ! Gustavson's algorithm for matrix-matrix multiplication
    ! THRESHOLDON: 0 = do not throw away elements; 1 = throw away elements
    ! NUMTHRESH: If THRESHOLDON = 1 throw away element whose absolute value is
    ! smaller than NUMTHRESH
    ! FILLINSTOP: Number of purification cycles beyond which we stop allowing
    ! for further fill-in
    !

    SPARSEON = VALVECTOR_INT(17); THRESHOLDON = VALVECTOR_INT(18); NUMTHRESH = VALVECTOR_RE(13)

#ifdef PROGRESSON

    IF(SPARSEON == 0)THEN
       LT%BML_TYPE = BML_MATRIX_DENSE
    ELSEIF(SPARSEON == 1)THEN
       LT%BML_TYPE = BML_MATRIX_ELLPACK
    ELSE
       CALL ERRORS("latteparser_latte_mod","SPARSEON > 1 yet not implemented")
    ENDIF

    IF(THRESHOLDON == 0)THEN
       LT%THRESHOLD = 0.0_DP
    ELSEIF(THRESHOLDON == 1)THEN
       LT%THRESHOLD = NUMTHRESH
    ELSE
       CALL ERRORS("latteparser_latte_mod","THRESHOLDON > 1 yet not implemented")
    ENDIF

#endif

    FILLINSTOP = VALVECTOR_INT(19); BLKSZ = VALVECTOR_INT(20)

    !
    ! MSPARSE: value for M when SPARSEON = 1, used by sp2 sparse algorithm
    !          0 = value for M is not known, defaults to N
    !

    MSPARSE = VALVECTOR_INT(21)

#ifdef PROGRESSON

    IF(MSPARSE == 0)THEN
       LT%MDIM = -1  !Defaults to N
    ELSEIF(MSPARSE > 0)THEN
       LT%MDIM = MSPARSE
    ELSE
       CALL ERRORS("latteparser_latte_mod","MSPARSE cannot be negative")
    ENDIF

#endif

    !
    ! LCNON: 0 = during charge neutral MD simulations we'll run LCNITER SCF
    ! cycles at each time step, 1 = we'll run SCF cycles until CHTOL is reached
    ! LCNITER: Number of SCF cycles to achieve LCN at each MD time step
    ! CHTOL: Tolerance on atomic charges (Mulliken) before LCN is declared
    !

    LCNON = VALVECTOR_INT(22); LCNITER = VALVECTOR_INT(23); CHTOL = VALVECTOR_RE(14)

    !
    ! Read the SKIN for the neighbor list (Angstrom)
    !

    SKIN = VALVECTOR_RE(15)

    !
    ! RELAXME: 0 = Don't run relaxation, 1 = relax geometry
    ! RELTYPE: SD = steepest descent, CG = conjugate gradient
    ! MXRLX: Maximum number of steps in the geometry optimization
    ! RLXFTOT: Run optimization until all forces are less than RLXFTOL
    !

    RELAXME = VALVECTOR_INT(24); RELTYPE = VALVECTOR_CHAR(4) ;MXRLX = VALVECTOR_INT(25)
    RLXFTOL = VALVECTOR_RE(16)

    !
    ! MDON: 0 = Molecular dynamics off, 1 = Molecular dynamics on
    ! (MD is controlled using the file MDcontroller)
    !

    MDON = VALVECTOR_INT(26)

    !
    ! PBCON: 1 = full periodic boundary conditions, 0 = gas phase: no pbc and
    ! electrostatics done all in real space
    !

    PBCON = VALVECTOR_INT(27)

    RESTART = VALVECTOR_INT(28)

    ! Add or remove electrons. 2+ -> charge = +2 since TOTNE = TOTNE - CHARGE

    CHARGE = VALVECTOR_INT(29)

    !
    ! XBOON: 0 = Niklasson's extended Lagrangian Born-Oppenheimer MD off,
    ! 1 = on.
    !

    XBOON = VALVECTOR_INT(30)

    !
    ! XBODISON: We have the option of turning on damping for the XBO
    ! to remedy the accumulation of noise. 0 = off, 1 = on.
    !

    XBODISON = VALVECTOR_INT(31)

    !
    ! XBODISORDER: = Order of the damping function (1 - 9)
    !

    XBODISORDER = VALVECTOR_INT(32)

    FITON = 0

    !
    ! Read in the number of GPUs per node
    !

    NGPU = VALVECTOR_INT(33)

    ! Are we doing k-space?

    KON = VALVECTOR_INT(34)

    ! Do we want to calculate forces too (not always necessary when fitting)

    COMPFORCE = VALVECTOR_INT(35)

    ! Turn on the simulated annealing subroutine to fit DOS


    DOSFITON = VALVECTOR_INT(36); INT2FIT = VALVECTOR_INT(37)
    MCBETA = VALVECTOR_RE(17); NFITSTEP = VALVECTOR_INT(38); QFIT = VALVECTOR_INT(39)
    MCSIGMA = VALVECTOR_RE(18)

    PPFITON = VALVECTOR_INT(40)

    ALLFITON = VALVECTOR_INT(41)

    PPNFITSTEP = VALVECTOR_INT(42); BINFITSTEP = VALVECTOR_INT(43)
    PP2FIT = VALVECTOR_INT(44); BINT2FIT = VALVECTOR_INT(45)

    PPBETA = VALVECTOR_RE(19); PPSIGMA = VALVECTOR_RE(20)
    PPNMOL = VALVECTOR_INT(46); PPNGEOM = VALVECTOR_INT(47)

    PARREP = VALVECTOR_INT(48)

    ! Verbosity level to control general output

    VERBOSE = VALVECTOR_INT(49)

    ! If Pulay Mixer

    MIXER = VALVECTOR_INT(50)

    ! Restart option for latte lib.

    MIXER = VALVECTOR_INT(51)

    ! Dielectric constant

    RELPERM = VALVECTOR_RE(21)

    ! Coordinates and parameter paths

    PARAMPATH = VALVECTOR_CHAR(5)
    COORDSFILE = VALVECTOR_CHAR(6)

    ! If latte_lib has to be restarted

    RESTARTLIB = VALVECTOR_INT(51)

    ! Freeze atoms

    FREEZE = VALVECTOR_INT(52)

    IF (ELECTRO == 1 .AND. ELECMETH == 1 .AND. PBCON == 1) THEN
      CALL ERRORS("latteparser_latte_mod","If CONTROL{ELECTRO= 1 ELECMETH= 1} &
      &then CONTROL{PBCON= 0}")
    END IF

    SCLTYPE = VALVECTOR_CHAR(7)

    STOPATMAXSCF = VALVECTOR_LOG(2)

  END SUBROUTINE PARSE_CONTROL


  !> The parser for Latte General input variables.
  !!
  SUBROUTINE PARSE_MD(FILENAME)

    IMPLICIT NONE
    INTEGER, PARAMETER :: NKEY_CHAR = 3, NKEY_INT = 18, NKEY_RE = 10, NKEY_LOG = 1
    CHARACTER(LEN=*) :: FILENAME

    !Library of keywords with the respective defaults.
    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'RNDIST=','SEEDINIT=','NPTTYPE=']
    CHARACTER(LEN=100) :: VALVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'GAUSSIAN','UNIFORM','ISO']

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_INT(NKEY_INT) = [CHARACTER(LEN=50) :: &
         'MAXITER=', 'UDNEIGH=', 'DUMPFREQ=','RSFREQ=', 'WRTFREQ=', 'TOINITTEMP5=', 'THERMPER=',& !7
         'THERMRUN=', 'NVTON=', 'NPTON=', 'AVEPER=', 'SEED=', 'SHOCKON=',&
         'SHOCKSTART=','SHOCKDIR=','MDADAPT=','GETHUG=','RSLEVEL=']
    INTEGER :: VALVECTOR_INT(NKEY_INT) = (/ &
         5000,1,250,500,25,1,500, &
         50000,0,0,1000,54,0, &
         100000,1,0,0,0/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_RE(NKEY_RE) = [CHARACTER(LEN=50) :: &
         'DT=','TEMPERATURE=','FRICTION=','PTARGET=','UPARTICLE=','USHOCK=','C0=', 'E0=',&
         'V0=','P0=']
    REAL(DP) :: VALVECTOR_RE(NKEY_RE) = (/&
         0.25,300.00,1000.0,0.0,500.0,-4590.0,1300.0,-795.725,&
         896.984864,0.083149/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_LOG(NKEY_LOG) = [CHARACTER(LEN=100) :: &
         'DUMMY=']
    LOGICAL :: VALVECTOR_LOG(NKEY_LOG) = (/&
         .FALSE./)

    !Start and stop characters
    CHARACTER(LEN=50), PARAMETER :: STARTSTOP(2) = [CHARACTER(LEN=50) :: &
         'MDCONTROL{', '}']

    CALL PARSING_KERNEL(KEYVECTOR_CHAR,VALVECTOR_CHAR&
         ,KEYVECTOR_INT,VALVECTOR_INT,KEYVECTOR_RE,VALVECTOR_RE,&
         KEYVECTOR_LOG,VALVECTOR_LOG,TRIM(FILENAME),STARTSTOP)


    !
    ! MAXITER = run this many MD time steps
    !

    MAXITER = VALVECTOR_INT(1)

    !
    ! UDNEIGH = update the neighbor lists every UDNEIGH time steps
    !

    UDNEIGH = VALVECTOR_INT(2)
    !   write(*,*)"UDNEIGH",UDNEIGH
    !
    ! DT = size of the time step in fs
    !

    DT = VALVECTOR_RE(1)
    !   write(*,*)"DT",DT
    !
    ! TTARGET = temperature in K were initialize and aim for during NVT MD
    ! RNDIST = Type of distribution of random numbers used to initialize T
    !        = GAUSSIAN or UNIFORM
    ! SEEDINIT = Type of seed used in the generation of random numbers
    !          = RANDOM - seed changes every time
    !          = DEFAULT - use the same, default seed every time
    !

    TTARGET = VALVECTOR_RE(2); RNDIST = VALVECTOR_CHAR(1); SEEDINIT = VALVECTOR_CHAR(2)
    !   write(*,*)"TTARGET,RNDIST,SEEDINIT",TTARGET,RNDIST,SEEDINIT
    !
    ! DUMPFREQ: Write a dump file every DUMPFREQ time steps
    !

    DUMPFREQ = VALVECTOR_INT(3)
    !   write(*,*)"DUMPFREQ",DUMPFREQ
    !
    ! RSFREQ: Write a restart file every RSFREQ time steps
    !

    RSFREQ = VALVECTOR_INT(4)
    !   write(*,*)"RSFREQ",RSFREQ
    !
    ! WRTFREQ: Output energy and temperature every WRTFREQ time steps
    !

    WRTFREQ = VALVECTOR_INT(5)
    IF ( WRTFREQ <= 0 ) THEN
      CALL ERRORS("latteparser_latte_mod","You cannot have WRTFREQ <= 0.&
                   &Set this variable to a very high value to avoid frequent printing")
    ENDIF
    !   write(*,*)"WRTFREQ",WRTFREQ
    !
    ! TOINITTEMP: Whether or not we are going to initialize velocities
    ! using a random number generator (sometimes during a restart we
    ! may not want to reinitialize the temperature
    !

    TOINITTEMP = VALVECTOR_INT(6)
    !   write(*,*)"TOINITTEMP",TOINITTEMP
    !
    ! THERMPER: If we're running NVT, rescale velocities every THERMPER
    ! time steps.
    !

    THERMPER = VALVECTOR_INT(7)
    !   write(*,*)"THERMPER",THERMPER
    !
    ! THERMRUN: Thermalize over this many time steps when NVT is on
    !

    THERMRUN = VALVECTOR_INT(8)
    !   write(*,*)"THERMRUN",THERMRUN
    !
    ! NVTON: 0 = running NVE MD, 1 = running NVT MD
    ! AVEPER: Average the temperature over AVEPER time steps when determining
    ! how to rescale velocities
    !

    NVTON = VALVECTOR_INT(9); NPTON = VALVECTOR_INT(10)
    !   write(*,*)"NVTON,NPTON",NVTON,NPTON
    AVEPER = VALVECTOR_INT(11); FRICTION = VALVECTOR_RE(3); SEEDTH = VALVECTOR_INT(12)
    !   write(*,*)"AVEPER,FRICTION,SEEDTH",AVEPER,FRICTION,SEEDTH


    IF (NVTON .EQ. 1 .AND. NPTON .EQ. 1) THEN
       CALL ERRORS("latteparser_latte_mod","You can't have NVTON = 1 and NPTON = 1")
    ENDIF

    ! PTARGET = Target pressure (in GPa) when running NPT
    ! NPTTYPE = ISO or ANISO

    PTARGET = VALVECTOR_RE(4); NPTTYPE = VALVECTOR_CHAR(3)
    !   write(*,*)"PTARGET,NPTTYPE",PTARGET,NPTTYPE

    !
    ! The following are for the Hugoniostat
    !

    ! On (1) or off (0)?

    SHOCKON = VALVECTOR_INT(13)
    !   write(*,*)"SHOCKON",SHOCKON
    !
    ! SHOCKSTART = the MD iteration where we will start to compress
    ! the iteration when we stop depends on the size of the block and Us
    !

    SHOCKSTART = VALVECTOR_INT(14)
    !   write(*,*)"SHOCKSTART",SHOCKSTART
    !
    ! SHOCKDIR is the cartensian direction (1 = X, 2 = Y, 3 = Z),
    ! parallel to which we're going to compress uniaxially
    !

    SHOCKDIR = VALVECTOR_INT(15)
    !   write(*,*)"SHOCKDIR",SHOCKDIR
    !
    ! And finally, the particle and shock velocities
    ! IN UNITS OF METRES PER SECOND
    !

    UPARTICLE = VALVECTOR_RE(5); USHOCK = VALVECTOR_RE(6); C0 = VALVECTOR_RE(7)
    !   write(*,*)"UPARTICLE,USHOCK,C0",UPARTICLE,USHOCK,C0
    ! Adapt SCF on the fly?

    MDADAPT = VALVECTOR_INT(16)
    !   write(*,*)"MDADAPT",MDADAPT
    ! Calculating Hugoniot points?

    GETHUG = VALVECTOR_INT(17)

    RSLEVEL = VALVECTOR_INT(18)

    !   write(*,*)"GETHUG",GETHUG
    E0 = VALVECTOR_RE(8); V0 = VALVECTOR_RE(9); P0 = VALVECTOR_RE(10)
    !   write(*,*)"E0,V0,P0",E0,V0,P0

  END SUBROUTINE PARSE_MD


  !> The parser for K Mesh input variables.
  !!
  SUBROUTINE PARSE_KMESH(FILENAME)

    IMPLICIT NONE
    INTEGER, PARAMETER :: NKEY_CHAR = 1, NKEY_INT = 3, NKEY_RE = 3, NKEY_LOG = 1
    CHARACTER(LEN=*) :: FILENAME

    !Library of keywords with the respective defaults.
    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'DUMMY=']
    CHARACTER(LEN=100) :: VALVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'DUMMY']

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_INT(NKEY_INT) = [CHARACTER(LEN=50) :: &
         'NKX=','NKY=','NKZ=']
    INTEGER :: VALVECTOR_INT(NKEY_INT) = (/ &
         1,1,1/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_RE(NKEY_RE) = [CHARACTER(LEN=50) :: &
         'KSHIFTX=','KSHIFTY=','KSHIFTZ=']
    REAL(DP) :: VALVECTOR_RE(NKEY_RE) = (/&
         0.0,0.0,0.0/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_LOG(NKEY_LOG) = [CHARACTER(LEN=100) :: &
         'DUMMY=']
    LOGICAL :: VALVECTOR_LOG(NKEY_LOG) = (/&
         .FALSE./)

    !Start and stop characters
    CHARACTER(LEN=50), PARAMETER :: STARTSTOP(2) = [CHARACTER(LEN=50) :: &
         'KMESH{', '}']

    CALL PARSING_KERNEL(KEYVECTOR_CHAR,VALVECTOR_CHAR&
         ,KEYVECTOR_INT,VALVECTOR_INT,KEYVECTOR_RE,VALVECTOR_RE,&
         KEYVECTOR_LOG,VALVECTOR_LOG,TRIM(FILENAME),STARTSTOP)

    NKX= VALVECTOR_INT(1); NKY= VALVECTOR_INT(2); NKZ=VALVECTOR_INT(3)
    KSHIFT(1)= VALVECTOR_RE(1); KSHIFT(2)= VALVECTOR_RE(2); KSHIFT(3)= VALVECTOR_RE(3)


  END SUBROUTINE PARSE_KMESH


END MODULE LATTEPARSER_LATTE_MOD
