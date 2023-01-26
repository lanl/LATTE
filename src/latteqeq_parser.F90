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


  !> The parser for Latte General input variables.
  !!
  SUBROUTINE PARSE_CONTROL(FILENAME)
    use qeq_parameters
    use KERNELPARSER_MOD
    IMPLICIT NONE
    INTEGER, PARAMETER :: NKEY_CHAR = 1, NKEY_INT = 1, NKEY_RE = 1, NKEY_LOG = 3
    CHARACTER(LEN=*) :: FILENAME

    !Library of keywords with the respective defaults.
    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'OUTFILE=']
    CHARACTER(LEN=100) :: VALVECTOR_CHAR(NKEY_CHAR) = [CHARACTER(LEN=100) :: &
         'log.latte']

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_INT(NKEY_INT) = [CHARACTER(LEN=50) :: &
         'VERBOSE=']
    INTEGER :: VALVECTOR_INT(NKEY_INT) = (/ &
         1 /)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_RE(NKEY_RE) = [CHARACTER(LEN=50) :: &
         'KBT='] !1
    REAL(8) :: VALVECTOR_RE(NKEY_RE) = (/&
         0.1/)

    CHARACTER(LEN=50), PARAMETER :: KEYVECTOR_LOG(NKEY_LOG) = [CHARACTER(LEN=100) :: &
         'EXACT=','PRINTCHARGES=','PRINTPOT=']
    LOGICAL :: VALVECTOR_LOG(NKEY_LOG) = (/&
         .FALSE.,.FALSE.,.FALSE./)

    !Start and stop characters
    CHARACTER(LEN=50), PARAMETER :: STARTSTOP(2) = [CHARACTER(LEN=50) :: &
         'CONTROL{', '}']

    CALL PARSING_KERNEL(KEYVECTOR_CHAR,VALVECTOR_CHAR&
         ,KEYVECTOR_INT,VALVECTOR_INT,KEYVECTOR_RE,VALVECTOR_RE,&
         KEYVECTOR_LOG,VALVECTOR_LOG,TRIM(FILENAME),STARTSTOP)

    VERBOSE = VALVECTOR_INT(1)

    exact_solution = VALVECTOR_LOG(1)
    printcharges = VALVECTOR_LOG(2)
    printpot = VALVECTOR_LOG(3)

  END SUBROUTINE PARSE_CONTROL

