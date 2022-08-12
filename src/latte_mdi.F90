MODULE POSIX
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_INT32_T
  IMPLICIT NONE

  INTERFACE
    ! INT USLEEP(USECONDS_T USECONDS)
    FUNCTION C_USLEEP(USECONDS) BIND(C, NAME='usleep')
      IMPORT :: C_INT, C_INT32_T
      INTEGER(KIND=C_INT32_T), VALUE :: USECONDS
      INTEGER(KIND=C_INT)            :: C_USLEEP
    END FUNCTION C_USLEEP
  END INTERFACE
END MODULE POSIX

! ------------------------------------------------------------------

MODULE LATTE_MDI

  USE LATTE_LIB

  USE MPI
  USE ISO_C_BINDING
  USE POSIX
  USE MDI, ONLY : MDI_INIT, MDI_SEND, MDI_INT, MDI_CHAR, MDI_NAME_LENGTH, &
       MDI_ACCEPT_COMMUNICATOR, MDI_RECV_COMMAND, MDI_RECV, &
       MDI_SET_EXECUTE_COMMAND_FUNC, MDI_MPI_GET_WORLD_COMM, MDI_DOUBLE, MDI_BYTE, &
       MDI_ENGINE, MDI_GET_ROLE, MDI_REGISTER_COMMAND, MDI_REGISTER_NODE, &
       MDI_REGISTER_CALLBACK, MDI_COMMAND_LENGTH, MDI_MPI_GET_WORLD_COMM, &
       MDI_PLUGIN_GET_ARG, MDI_PLUGIN_GET_ARGC

  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)

  ! MDI Communicator to the driver
  INTEGER :: COMM

  ! MPI intra-communicator for this code
  INTEGER :: WORLD
  INTEGER :: ME, NPROCS

  ! Flag to terminate MDI response function
  LOGICAL :: EXITFLAG = .FALSE.

  ! exchanged data
  INTEGER :: NATOMS, NTYPES, LFNAME
  INTEGER, ALLOCATABLE :: TYPES(:)
  CHARACTER(100) :: SYMBLIST
  CHARACTER(100) :: FNAME = "latte.in"
  DOUBLE PRECISION, ALLOCATABLE :: AUX(:), CELL(:), CELL_DISPL(:)
  DOUBLE PRECISION :: TDELAY, VENERG

  ! Data that needs to stay in memory
  CHARACTER(2), ALLOCATABLE :: SYMB(:)
  INTEGER, ALLOCATABLE :: ELEMENTS(:)
  DOUBLE PRECISION, ALLOCATABLE :: COORDS(:,:), MASSES(:), FORCES(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: STRESS(:), BOX(:,:)
  INTEGER :: NEWSYSTEM = 0
  INTEGER :: MDIVERB = 0
  CHARACTER(1) :: STATE
  LOGICAL :: FIRSTSYSTEM = .true.

  !> Element symbols
  !!
  character(2), parameter :: MDI_element_symbol(103) = [character(2) :: &
       "H" ,          "He" ,         "Li" ,         "Be" ,         &
       "B" ,          "C" ,          "N" ,          "O" ,          &
       "F" ,          "Ne" ,         "Na" ,         "Mg" ,         &
       "Al" ,         "Si" ,         "P" ,          "S" ,          &
       "Cl" ,         "Ar" ,         "K" ,          "Ca" ,         &
       "Sc" ,         "Ti" ,         "V" ,          "Cr" ,         &
       "Mn" ,         "Fe" ,         "Co" ,         "Ni" ,         &
       "Cu" ,         "Zn" ,         "Ga" ,         "Ge" ,         &
       "As" ,         "Se" ,         "Br" ,         "Kr" ,         &
       "Rb" ,         "Sr" ,         "Y" ,          "Zr" ,         &
       "Nb" ,         "Mo" ,         "Tc" ,         "Ru" ,         &
       "Rh" ,         "Pd" ,         "Ag" ,         "Cd" ,         &
       "In" ,         "Sn" ,         "Sb" ,         "Te" ,         &
       "I" ,          "Xe" ,         "Cs" ,         "Ba" ,         &
       "La" ,         "Ce" ,         "Pr" ,         "Nd" ,         &
       "Pm" ,         "Sm" ,         "Eu" ,         "Gd" ,         &
       "Tb" ,         "Dy" ,         "Ho" ,         "Er" ,         &
       "Tm" ,         "Yb" ,         "Lu" ,         "Hf" ,         &
       "Ta" ,         "W" ,          "Re" ,         "Os" ,         &
       "Ir" ,         "Pt" ,         "Au" ,         "Hg" ,         &
       "Tl" ,         "Pb" ,         "Bi" ,         "Po" ,         &
       "At" ,         "Rn" ,         "Fr" ,         "Ra" ,         &
       "Ac" ,         "Th" ,         "Pa" ,         "U" ,          &
       "Np" ,         "Pu" ,         "Am" ,         "Cm" ,         &
       "Bk" ,         "Cf" ,         "Es" ,         "Fm" ,         &
       "Md" ,         "No" ,         "Lr"                          &
       ]

CONTAINS

  SUBROUTINE GET_TYPES_AND_SYMBOLS(ELEMENTS,TYPES,SYMBOLS)
    IMPLICIT NONE 
    INTEGER, ALLOCATABLE, INTENT(IN) :: ELEMENTS(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: TYPES(:)
    CHARACTER(2), ALLOCATABLE, INTENT(OUT) :: SYMBOLS(:)
    integer :: nats, tmpcount,i,j
    integer, allocatable :: tmptypes(:)

    NATS = SIZE(ELEMENTS)
    ALLOCATE(TMPTYPES(NATS)) !Temporary allocation
    TMPTYPES = 0

    !Count the types
    NTYPES = 1
    TMPTYPES(1) = ELEMENTS(1)
    DO I = 1,NATS 
      TMPCOUNT = 0
      DO J = 1,NATS 
        IF(TMPTYPES(J) == 0) EXIT
        IF(ELEMENTS(I) /= TMPTYPES(J))THEN
          TMPCOUNT = TMPCOUNT + 1
        ENDIF
      ENDDO
      IF(TMPCOUNT >= NTYPES)THEN !There is a new type
        NTYPES = NTYPES + 1
        TMPTYPES(NTYPES) = ELEMENTS(I)
      ENDIF
    ENDDO

    ALLOCATE(TYPES(NATS))
    do i = 1,NATS
      do j = 1,NTYPES
        if(ELEMENTS(i) == TMPTYPES(j))THEN
          TYPES(i) = j
        endif
      enddo
    enddo
    ALLOCATE(SYMBOLS(NTYPES))
    DO I = 1,NTYPES
     SYMBOLS(I) = MDI_ELEMENT_SYMBOL(TMPTYPES(I)) 
    ENDDO

    DEALLOCATE(TMPTYPES)

  END SUBROUTINE GET_TYPES_AND_SYMBOLS
  ! -----------------------------------------------------------------

  FUNCTION MDI_PLUGIN_INIT_LATTE_MDI() BIND (C, NAME="MDI_Plugin_init_latte_mdi")
    INTEGER :: MDI_PLUGIN_INIT_LATTE_MDI
    INTEGER :: IERR
    INTEGER :: ARGC
    INTEGER :: IARG
    CHARACTER(LEN=1024) :: OPTION
    CHARACTER(LEN=1024) :: MDI_OPTION, OTHER_OPTIONS
    LOGICAL :: MDI_OPTION_FOUND
    
    ! how to declare other_options as vector of strings?
    MDI_OPTION_FOUND = .FALSE.
    CALL MDI_PLUGIN_GET_ARGC(ARGC,IERR)

    DO IARG=0, ARGC-1
      CALL MDI_PLUGIN_GET_ARG(IARG,OPTION,IERR)
      IF ( (TRIM(OPTION) .EQ. "-mdi") .OR. (TRIM(OPTION) .EQ. "--mdi") ) THEN
        IF ( ARGC .GT. (IARG+1) ) THEN
          CALL MDI_PLUGIN_GET_ARG(IARG+1, MDI_OPTION, IERR)
          MDI_OPTION_FOUND = .TRUE.
        ELSE
          WRITE(6,*)'ERROR: LATTE -mdi argument not provided'
          MDI_PLUGIN_INIT_LATTE_MDI = 1
          RETURN
        END IF
      ELSE
        ! how to copy arg into other_options vector
      END IF
    END DO

    IF (.NOT. MDI_OPTION_FOUND) THEN
      WRITE(6,*)'ERROR: LATTE -mdi option not provided'
      MDI_PLUGIN_INIT_LATTE_MDI = 1
      RETURN
    END IF

    ! start LATTE running as an MDI engine
    CALL MDI_ENGINE_INVOKE(MDI_OPTION,OTHER_OPTIONS)

    MDI_PLUGIN_INIT_LATTE_MDI = 0

  END FUNCTION MDI_PLUGIN_INIT_LATTE_MDI

  ! -----------------------------------------------------------------

  SUBROUTINE OPTIONS(OTHER_OPTIONS)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: OTHER_OPTIONS
    INTEGER :: IERR

    ! print out list of other option flags

    WRITE (6,*) "OPTIONS"

  END SUBROUTINE OPTIONS

  ! -----------------------------------------------------------------

  SUBROUTINE MDI_ENGINE_INVOKE(MDI_OPTION,OTHER_OPTIONS)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: MDI_OPTION, OTHER_OPTIONS
    INTEGER :: IERR, ROLE
    CHARACTER(LEN=:), ALLOCATABLE :: COMMAND

    PROCEDURE(EXECUTE_COMMAND), POINTER :: GENERIC_COMMAND => NULL()
    TYPE(C_PTR)                         :: CLASS_OBJ
    GENERIC_COMMAND => EXECUTE_COMMAND

    ! Call MDI_Init

    CALL MDI_INIT(MDI_OPTION,IERR)

    ! Get the MPI intra-communicator over which this plugin will run

    CALL MDI_MPI_GET_WORLD_COMM(WORLD,IERR)
    CALL MPI_COMM_RANK(WORLD, ME, IERR)
    CALL MPI_COMM_SIZE(WORLD, NPROCS, IERR)

    ! process non-MDI command line args

    CALL OPTIONS(OTHER_OPTIONS)

    ! Confirm LATTE is being run as an ENGINE

    CALL MDI_GET_ROLE(ROLE, IERR)
    IF ( ROLE .NE. MDI_ENGINE ) THEN
      WRITE(6,*)'ERROR: Must run engine_f90 as an ENGINE'
    END IF

    ! supported MDI commands

    CALL MDI_REGISTER_NODE("@DEFAULT", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "EXIT", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">FNAME", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">VERBOSE", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">NATOMS", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">ELEMENTS", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">COORDS", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">CELL", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", ">CELL_DISPL", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "RUN", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "NEWSYSTEM", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "<FORCES", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "<PE", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "<STRESS", IERR)
    CALL MDI_REGISTER_COMMAND("@DEFAULT", "ECHO", IERR)
    
    ! one-time operation to establish a connection WITH the driver

    CALL MDI_ACCEPT_COMMUNICATOR(COMM, IERR)

    ! set callback to execute_command

    CALL MDI_SET_EXECUTE_COMMAND_FUNC(GENERIC_COMMAND, CLASS_OBJ, IERR)

    ALLOCATE( CHARACTER(MDI_COMMAND_LENGTH) :: COMMAND )

    ! Respond to the driver's commands

    RESPONSE_LOOP: DO

      ! Receive a command from the driver and broadcast it to all ranks
      CALL MDI_RECV_COMMAND(COMMAND, COMM, IERR)
      CALL MPI_BCAST(COMMAND, MDI_COMMAND_LENGTH, MPI_CHAR, 0, WORLD, IERR)
      CALL EXECUTE_COMMAND(COMMAND, COMM, IERR)
      IF (EXITFLAG) EXIT

    END DO RESPONSE_LOOP

    IF(ALLOCATED(TYPES)) DEALLOCATE(TYPES)
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
    IF(ALLOCATED(CELL)) DEALLOCATE(CELL)
    IF(ALLOCATED(SYMB)) DEALLOCATE(SYMB)
    IF(ALLOCATED(ELEMENTS)) DEALLOCATE(ELEMENTS)
    IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
    IF(ALLOCATED(MASSES)) DEALLOCATE(MASSES)
    IF(ALLOCATED(FORCES)) DEALLOCATE(FORCES)
    IF(ALLOCATED(COMMAND)) DEALLOCATE(COMMAND)

  END SUBROUTINE MDI_ENGINE_INVOKE

  ! -----------------------------------------------------------------

  SUBROUTINE EXECUTE_COMMAND(COMMAND, MDICOMM, IERR)
    USE LATTE_LIB
    USE SETUPARRAY, ONLY : NORBINDEX
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: COMMAND
    INTEGER, INTENT(IN)           :: MDICOMM
    INTEGER, INTENT(OUT)          :: IERR

    INTEGER :: US,RC

    !Latte_lib specifics
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER ::  MAXITER, I, J, IOS
    LOGICAL(1) :: EXISTERROR
    REAL(DP), ALLOCATABLE ::  VEL(:,:), CHARGES(:)
    CHARACTER(2), ALLOCATABLE :: SYMBOLS(:)
    REAL(DP) :: VENERGTMP
    REAL(dp), ALLOCATABLE :: STRESSTMP(:)

    INTEGER, PARAMETER :: LNAME = 100
    REAL(DP), PARAMETER :: eV2H = 1.0_dp/27.211396641308_dp
    REAL(dp), PARAMETER :: Ang2Bohr = 1.0_dp/0.529177_dp


    IF(MDIVERB >= 1)THEN 
      WRITE(*,*)"MDISTATE = ",STATE
      WRITE(*,*)"FIRSTSYSYTEM = ",FIRSTSYSTEM
      WRITE(*,*)"NEWSYSTEM = ",NEWSYSTEM
    ENDIF 

    SELECT CASE(TRIM(COMMAND))

    CASE("EXIT")
      EXITFLAG = .TRUE.

    ! Receving the name of the latte file
    CASE(">FNAME")
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving FNAME ..."
      FNAME = ""
      CALL MDI_RECV(FNAME, LNAME, MDI_CHAR, MDICOMM, IERR)
      IF(IERR > 0)THEN
        WRITE(*,*)"ERROR passing FNAME. Check the leght name integer&
             & passed from the driver. The length integer should be set to: ",LNAME
        STOP
      ENDIF
      CALL MPI_BCAST(FNAME, LNAME, MPI_CHAR, 0, WORLD, IERR)
      !WRITE(*,*)"Name of latte file ",FNAME

    ! Receving verbosity level 
    CASE(">VERBOSE")
      WRITE(*,*)"Receiving VERBOSE ..."
      WRITE(*,*)"VERBOSE is ON ..."
      CALL MDI_RECV(MDIVERB, 1, MDI_INT, MDICOMM, IERR)
      CALL MPI_BCAST(MDIVERB, 1, MPI_INT, 0, WORLD, IERR)
      IF(IERR > 0)THEN
        WRITE(*,*)"ERROR passing FNAME. Check the leght name integer&
             & passed from the driver. The length integer should be set to: ",LNAME
        STOP
      ENDIF
      CALL MPI_BCAST(FNAME, LNAME, MPI_CHAR, 0, WORLD, IERR)
      !WRITE(*,*)"Name of latte file ",FNAME
  

    ! Receiving the number of atoms
    CASE( ">NATOMS" )
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving NATOMS ..."
      CALL MDI_RECV(NATOMS, 1, MDI_INT, MDICOMM, IERR)
      CALL MPI_BCAST(NATOMS, 1, MPI_INT, 0, WORLD, IERR)
      IF(ALLOCATED(ELEMENTS)) DEALLOCATE(ELEMENTS)
      IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
      IF(ALLOCATED(FORCES)) DEALLOCATE(FORCES)
      IF(ALLOCATED(NORBINDEX)) DEALLOCATE(NORBINDEX)
      IF(MDIVERB >= 1) WRITE(*,*)"Received ",NATOMS," atoms ..."
      IF(FIRSTSYSTEM)THEN 
        NEWSYSTEM = 0
        FIRSTSYSTEM = .false.
      ELSE
        NEWSYSTEM = 1
      ENDIF
      STATE = ">"

    ! Receiving element atomic numbers
    CASE( ">ELEMENTS" )
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving ELEMENTS ..."
      IF(NEWSYSTEM == 1)THEN 
        IF(ALLOCATED(ELEMENTS))DEALLOCATE(ELEMENTS)
      ENDIF
      IF(.NOT. ALLOCATED(ELEMENTS)) ALLOCATE(ELEMENTS(NATOMS))
      CALL MDI_RECV(ELEMENTS, NATOMS, MDI_INT, MDICOMM, IERR)
      CALL MPI_BCAST(ELEMENTS, NATOMS, MPI_INT, 0, WORLD, IERR)
      CALL GET_TYPES_AND_SYMBOLS(ELEMENTS,TYPES,SYMB)
      IF(FIRSTSYSTEM)THEN
        NEWSYSTEM = 0
        FIRSTSYSTEM = .false.
      ELSE
        NEWSYSTEM = 1
      ENDIF
      STATE = ">"

    ! Receiving the coordinate. A 3*nats auxiliary array is used
    ! to pass the coordinated.
    CASE( ">COORDS" )
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving COORDS"
      ALLOCATE(AUX(3*NATOMS))
      IF(.NOT. ALLOCATED(COORDS)) ALLOCATE(COORDS(3,NATOMS))
      CALL MDI_RECV(AUX, 3*NATOMS, MDI_DOUBLE, MDICOMM, IERR)
      CALL MPI_BCAST(AUX, 3*NATOMS, MPI_DOUBLE, 0, WORLD, IERR)
      AUX = AUX/Ang2Bohr
      DO I = 1, NATOMS
        COORDS(1,I) = AUX(3*(I-1)+1)
        COORDS(2,I) = AUX(3*(I-1)+2)
        COORDS(3,I) = AUX(3*(I-1)+3)
      ENDDO
      DEALLOCATE(AUX)
      STATE = ">"

    ! Receiving the cell. The format that is passed is the same
    ! as the one used by lammps.
    CASE( ">CELL" )
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving CELL ..."
      ALLOCATE(CELL(9))
      IF(.not.ALLOCATED(BOX)) ALLOCATE(BOX(3,3))
      CALL MDI_RECV(CELL, 9, MDI_DOUBLE, MDICOMM, IERR)
      CALL MPI_BCAST(CELL, 9, MPI_DOUBLE, 0, WORLD, IERR)
      CELL = CELL/Ang2Bohr
      BOX(1,1) = CELL(1); BOX(1,2) = CELL(2); BOX(1,3) = CELL(3)
      BOX(2,1) = CELL(4); BOX(2,2) = CELL(5); BOX(2,3) = CELL(6)
      BOX(3,1) = CELL(7); BOX(3,2) = CELL(8); BOX(3,3) = CELL(9)
      DEALLOCATE(CELL)
      STATE = ">"

    ! Receiving the cell displacement.
    CASE( ">CELL_DISPL" )
      IF(MDIVERB >= 1) WRITE(*,*)"Receiving CELL_DISP ..."
      ALLOCATE(CELL_DISPL(3))
      CELL_DISPL = CELL_DISPL/Ang2Bohr
      CALL MDI_RECV(CELL_DISPL, 3, MDI_DOUBLE, MDICOMM, IERR)
      CALL MPI_BCAST(CELL_DISPL, 3, MPI_DOUBLE, 0, WORLD, IERR)
      IF(MDIVERB >= 1) WRITE(*,*)"WARNING: CELL_DISL is not used within LATTE"
      DEALLOCATE(CELL_DISPL)
      STATE = ">"

    ! This command will run latte, with the info received
    CASE( "RUN" )
      IF(MDIVERB >= 1) WRITE(*,*)"Executing RUN command ..."
      MAXITER = -1
      IF(.NOT. ALLOCATED(FORCES)) ALLOCATE(FORCES(3,NATOMS))
      FORCES = 0.0_DP
      iF(.NOT. ALLOCATED(STRESS)) ALLOCATE(STRESS(9))
      CALL LATTE(NTYPES, TYPES, COORDS, BOX, FORCES, &
           MAXITER, VENERG, STRESS, NEWSYSTEM, EXISTERROR, SYMB, FNAME)
      NEWSYSTEM = 0

    ! This command will clear the latte memory and recompute all the
    ! matrices
    CASE ( "NEWSYSTEM" )
      NEWSYSTEM = 1
      IF(MDIVERB >= 1) WRITE(*,*)"Setting up a newsystem ..."

    ! Passing back the forces using a 3*nats auxiliary array.
    CASE ( "<FORCES" )
      IF(MDIVERB >= 1) WRITE(*,*)"Sending FORCES"
      MAXITER = -1
      IF(.NOT. ALLOCATED(ELEMENTS)) STOP "ERROR: ELEMENTS were not received yet ..."
      IF(.NOT. ALLOCATED(COORDS)) STOP "ERROR: COORDS were not received yet ..."
      IF(.NOT. ALLOCATED(BOX)) STOP "ERROR: CELL was not received yet ..."
      IF(.NOT. ALLOCATED(FORCES)) ALLOCATE(FORCES(3,NATOMS))
      IF(.NOT. ALLOCATED(STRESS)) ALLOCATE(STRESS(9))
      IF ( STATE == ">" ) THEN
        FORCES = 0.0d0 
        STRESS = 0.0d0 
        VENERG = 0.0d0 
        IF(MDIVERB >= 2) WRITE(*,*)"BOX",BOX
        CALL LATTE(NTYPES, TYPES, COORDS, BOX, FORCES, &
        MAXITER, VENERG, STRESS, NEWSYSTEM, EXISTERROR, SYMB, FNAME)
        IF(MDIVERB >= 2) WRITE(*,*)"FORCESBACK",FORCES(1:3,1),"...",FORCES(1:3,NATOMS)
        NEWSYSTEM = 0
        ALLOCATE(AUX(3*NATOMS))
        DO I = 1,NATOMS
          AUX(3*(I-1) + 1) = FORCES(1,I)
          AUX(3*(I-1) + 2) = FORCES(2,I)
          AUX(3*(I-1) + 3) = FORCES(3,I)
        ENDDO
        AUX = AUX*eV2H/Ang2Bohr   !eV/Ang to H/Bohr
        IF(MDIVERB >= 2) WRITE(*,*)"FORCESAUX",AUX(1:3),"...",AUX(3*(NATOMS-1)+1:3*NATOMS)
        CALL MDI_SEND(AUX, 3*NATOMS, MDI_DOUBLE, MDICOMM, IERR)
        STATE = "<" 
      ELSE
        IF(MDIVERB >= 2) WRITE(*,*)"FORCESBACK",FORCES(1:3,1),"...",FORCES(1:3,NATOMS)
        ALLOCATE(AUX(3*NATOMS))
        DO I = 1,NATOMS
          AUX(3*(I-1) + 1) = FORCES(1,I)
          AUX(3*(I-1) + 2) = FORCES(2,I)
          AUX(3*(I-1) + 3) = FORCES(3,I)
        ENDDO
        AUX = AUX*eV2H/Ang2Bohr   !eV/Ang to H/Bohr
        IF(MDIVERB >= 2) WRITE(*,*)"FORCESAUX",AUX(1:3),"...",AUX(3*(NATOMS-1)+1:3*NATOMS) 
        CALL MDI_SEND(AUX, 3*NATOMS, MDI_DOUBLE, MDICOMM, IERR)
      ENDIF
      DEALLOCATE(AUX)

    ! Passing the potential energy
    CASE ( "<PE" )
      IF(MDIVERB >= 1) WRITE(*,*)"Sending PE ..."
      MAXITER = -1
      IF(.NOT. ALLOCATED(ELEMENTS)) STOP "EROOR: ELEMENTS were not received yet ..."
      IF(.NOT. ALLOCATED(COORDS)) STOP "EROOR: COORDS were not received yet ..."
      IF(.NOT. ALLOCATED(BOX)) STOP "EROOR: CELL was not received yet ..."
      IF(.NOT. ALLOCATED(FORCES)) ALLOCATE(FORCES(3,NATOMS))
      iF(.NOT. ALLOCATED(STRESS)) ALLOCATE(STRESS(9))
      IF ( STATE == ">" ) THEN
        FORCES = 0.0d0
        STRESS = 0.0d0
        VENERG = 0.0d0
        CALL LATTE(NTYPES, TYPES, COORDS, BOX, FORCES, &
        MAXITER, VENERG, STRESS, NEWSYSTEM, EXISTERROR, SYMB, FNAME)
        IF(MDIVERB >= 2) WRITE(*,*)"FORCESBACK",FORCES(1:3,1),"...",FORCES(1:3,NATOMS)
        NEWSYSTEM = 0
        VENERGTMP = VENERG*eV2H
        CALL MDI_SEND(VENERGTMP, 1, MDI_DOUBLE, MDICOMM, IERR)
        STATE = "<"
      ELSE
        VENERGTMP = VENERG*eV2H
        CALL MDI_SEND(VENERGTMP, 1, MDI_DOUBLE, MDICOMM, IERR)
      ENDIF


    ! Passing stress tensor 
    CASE ( "<STRESS" )
      MAXITER = -1
      IF(.NOT. ALLOCATED(ELEMENTS)) STOP "EROOR: ELEMENTS were not received yet ..."
      IF(.NOT. ALLOCATED(COORDS)) STOP "EROOR: COORDS were not received yet ..."
      IF(.NOT. ALLOCATED(BOX)) STOP "EROOR: CELL was not received yet ..."
      IF(.NOT. ALLOCATED(FORCES)) ALLOCATE(FORCES(3,NATOMS))
      iF(.NOT. ALLOCATED(STRESS)) ALLOCATE(STRESS(9))
      IF ( STATE == ">" ) THEN
        FORCES = 0.0d0
        STRESS = 0.0d0 
        VENERG = 0.0d0
        CALL LATTE(NTYPES, TYPES, COORDS, BOX, FORCES, &
        MAXITER, VENERG, STRESS, NEWSYSTEM, EXISTERROR, SYMB, FNAME)
        NEWSYSTEM = 0
        ALLOCATE(STRESSTMP(9))
        STRESSTMP = STRESS*eV2H/(Ang2Bohr**3)
        IF(MDIVERB >= 2) WRITE(*,*)"STRESS",STRESSTMP
        CALL MDI_SEND(STRESSTMP, 9, MDI_DOUBLE, MDICOMM, IERR)
        DEALLOCATE(STRESSTMP)  
        STATE = "<"
      ELSE
        ALLOCATE(STRESSTMP(9))
        STRESSTMP = STRESS*eV2H/(Ang2Bohr**3)
        IF(MDIVERB >= 2) WRITE(*,*)"STRESS",STRESSTMP
        CALL MDI_SEND(STRESSTMP, 9, MDI_DOUBLE, MDICOMM, IERR)
        DEALLOCATE(STRESSTMP)  
      ENDIF
      !write(*,*)"STRESS tensor", STRESS
      
    ! Print out the variables that were received
    CASE( "ECHO" )
      write(*,*)"Name of latte file ", FNAME
      write(*,*)"Number of atoms ", NATOMS
      write(*,*)"Number of species/atom types",NTYPES
      write(*,*)"List of elements in the system ",SYMB
      write(*,*)"List of Types in system ",TYPES
      write(*,*)"System coordinates ",COORDS
      write(*,*)"System cell xlo1, xhi1, xlo2, xhi2, xlo3, xhi3, xy, xz, yz ", CELL

    CASE DEFAULT
      WRITE(6,*)'ERROR: Unrecognized MDI command'

    END SELECT

    IERR = 0

  END SUBROUTINE EXECUTE_COMMAND

END MODULE LATTE_MDI

! ------------------------------------------------------------------
! main program
! ------------------------------------------------------------------

PROGRAM LATTE_MDI_ENGINE !Program

  USE CONSTANTS_MOD
  USE TIMER_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE SPARSEARRAY
  USE MDARRAY
  USE MYPRECISION
  USE VIRIALARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE LATTEPARSER
  USE NEBLISTARRAY
  USE NONOARRAY
  USE CONSTRAINTS_MOD

  USE LATTE_LIB

  USE MPI
  USE LATTE_MDI, ONLY : MDI_ENGINE_INVOKE

  IMPLICIT NONE

  INTEGER :: IARG, NARG, IERR
  CHARACTER(LEN=1024) :: ARG, MDI_OPTION, OTHER_OPTIONS
  LOGICAL :: MDI_OPTION_FOUND
  ! how to declare other_options as vector of strings?

  ! Initialize the MPI environment

  CALL MPI_INIT(IERR)

  ! mdi_option = single arg in quotes that follows -mdi
  ! other_options = all non-MDI args

  MDI_OPTION_FOUND = .FALSE.
  NARG = COMMAND_ARGUMENT_COUNT()

  IARG = 1
  DO WHILE (IARG <= NARG)
    CALL GET_COMMAND_ARGUMENT(IARG,ARG)

    IF (TRIM(ARG) .EQ. "-mdi" .OR. TRIM(ARG) == "--mdi") THEN
      IF (NARG >= IARG+1) THEN
        CALL GET_COMMAND_ARGUMENT(IARG+1,MDI_OPTION)
        MDI_OPTION_FOUND = .TRUE.
      ELSE
        WRITE(6,*) 'ERROR: LATTE -mdi argument not provided'
        STOP
      END IF
      IARG = IARG + 1
    ELSE
      ! how to copy arg into other_options vector
    END IF
    IARG = IARG + 1
  END DO

  IF (.NOT. MDI_OPTION_FOUND) THEN
    WRITE(6,*)'ERROR: LATTE -mdi option not provided'
    STOP
  END IF

  CALL MDI_ENGINE_INVOKE(MDI_OPTION,OTHER_OPTIONS)

  CALL MPI_FINALIZE(IERR)

END PROGRAM LATTE_MDI_ENGINE
