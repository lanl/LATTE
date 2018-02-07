!> Some general parsing functions.
!! \ingroup PROGRESS
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
MODULE KERNELPARSER_MOD

  USE OPENFILES_MOD
  USE CONSTANTS_MOD, ONLY: VERBOSE
  !  USE PARALLEL_MOD

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  PUBLIC :: PARSING_KERNEL

CONTAINS

  !> The general parsing function.
  !! It is used to vectorize a set of "keywords" "value" pairs
  !! as included in a general input file.
  !! \note This parsing strategy can only parse a file of
  !! 500 lines by 500 words.
  !! \warning If the length of variable vect is changed, this could produce a
  !! segmentation fault.
  !!
  SUBROUTINE PARSING_KERNEL(KEYVECTOR_CHAR,VALVECTOR_CHAR&
       ,KEYVECTOR_INT,VALVECTOR_INT,KEYVECTOR_RE,VALVECTOR_RE,&
       KEYVECTOR_LOG,VALVECTOR_LOG,FILENAME,STARTSTOP)
    IMPLICIT NONE
    CHARACTER(1), ALLOCATABLE              ::  tempc(:)
    CHARACTER(100), ALLOCATABLE            ::  vect(:,:)
    CHARACTER(50)                          ::  keyvector_char(:), keyvector_int(:), keyvector_log(:), keyvector_re(:)
    CHARACTER(100)                         ::  valvector_char(:)
    CHARACTER(LEN=*)                       ::  FILENAME
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL ::  STARTSTOP(2)
    CHARACTER(LEN=100)                      ::  TEMPCFLEX
    INTEGER                                ::  i, io_control, ios, j
    INTEGER                                ::  k, l, lenc, nkey_char
    INTEGER                                ::  nkey_int, nkey_log, nkey_re, readmaxi
    INTEGER                                ::  readmaxj, readmini, valvector_int(:)
    INTEGER                                ::  startatj, totalwords
    LOGICAL                                ::  start, stopl, valid, valvector_log(:), stopparsing, defaultnone
    LOGICAL, ALLOCATABLE                   ::  checkmissing_char(:), checkmissing_int(:), checkmissing_log(:), checkmissing_re(:)
    REAL(dp)                               ::  valvector_re(:)

    READMAXI = 5 ; READMAXJ = 1000
    ALLOCATE(VECT(READMAXI,READMAXJ))
    NKEY_CHAR = SIZE(KEYVECTOR_CHAR,DIM=1)
    NKEY_RE = SIZE(KEYVECTOR_RE,DIM=1)
    NKEY_INT = SIZE(KEYVECTOR_INT,DIM=1)
    NKEY_LOG = SIZE(KEYVECTOR_LOG,DIM=1)

    CALL OPEN_FILE_TO_READ(IO_CONTROL,FILENAME)

    ALLOCATE(CHECKMISSING_CHAR(NKEY_CHAR),CHECKMISSING_RE(NKEY_RE), &
         CHECKMISSING_INT(NKEY_INT), CHECKMISSING_LOG(NKEY_LOG))

    !Initialize the checkmissing flags and the vect arrays
    CHECKMISSING_CHAR = .FALSE.
    CHECKMISSING_RE = .FALSE.
    CHECKMISSING_INT = .FALSE.
    CHECKMISSING_LOG = .FALSE.
    STOPPARSING = .FALSE.
    DEFAULTNONE = .FALSE.
    VECT = '                    '

    DO I=1,READMAXI !Here we read all the input into vect
       READ(IO_CONTROL,*,IOSTAT=IOS)(VECT(I,J),J=1,READMAXJ)
    END DO

    CLOSE(IO_CONTROL)

    !Look up for floating hashes (#)
    TOTALWORDS = 0
    DO I=1,READMAXI
       DO K=1,READMAXJ
          IF(ADJUSTL(TRIM(VECT(I,K))).NE."")TOTALWORDS = TOTALWORDS + 1
          IF(ADJUSTL(TRIM(VECT(I,K))).EQ."#")THEN
             WRITE(*,*)" "
             WRITE(*,*)"ERROR in the the input file ..."
             WRITE(*,*)" "
             WRITE(*,*)"In the LATTE parsing routine everything is a comment by default unless theres an = sign"
             WRITE(*,*)"next to a word, in which case, it will be recognized as a keyword."
             WRITE(*,*)"This parser does not accept floating hashes (#). This is done in order to make sure"
             WRITE(*,*)"that a specific keyword is commented"
             WRITE(*,*)" "
             WRITE(*,*)"If you have a commented keyword make sure there is a # symbol right next to it"
             WRITE(*,*)"   "
             WRITE(*,*)"   The following commented keyword is correct: "
             WRITE(*,*)"                #KeyWord= 1 "
             WRITE(*,*)" "
             WRITE(*,*)"   The following commented keyword is NOT correct: "
             WRITE(*,*)"                # KeyWord= 1 "
             WRITE(*,*)" "
             STOP
          ENDIF
          IF(ADJUSTL(TRIM(VECT(I,K))).EQ."STOP{}")STOPPARSING = .TRUE.
          IF(ADJUSTL(TRIM(VECT(I,K))).EQ."DEFAULTNONE")DEFAULTNONE = .TRUE.
       ENDDO
    ENDDO

    IF(TOTALWORDS > READMAXI*READMAXJ - 100) THEN
       WRITE(*,*)""; WRITE(*,*)"Stopping ... Maximum allowed (keys + values + comments) words close to the limit "
       WRITE(*,*)"Increase the readmaxj variable in the parsing_kernel subroutine or reduce the comments in the input"
       STOP
    ENDIF

    !Look up for boundaries
    READMINI=1
    START=.FALSE.
    IF(PRESENT(STARTSTOP))THEN
       DO I=1,READMAXI
          DO K=1,READMAXJ
             IF(TRIM(VECT(I,K)).EQ.TRIM(STARTSTOP(1)))THEN
                READMINI=I
                STARTATJ=K
                START=.TRUE.
             ENDIF
             IF(START.AND.TRIM(VECT(I,K)).EQ.TRIM(STARTSTOP(2)))THEN
                READMAXI=I
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Look for invalid characters if startstop is present
    IF(START)THEN
       START=.FALSE.
       STOPL=.FALSE.
       DO I=READMINI,READMAXI
          DO K=1,READMAXJ
             IF(TRIM(VECT(I,K)).EQ.TRIM(STARTSTOP(1)))START=.TRUE.
             VALID = .FALSE.
             IF(START)THEN
                IF(VECT(I,K).NE.'                    ')THEN
                   DO J=1,NKEY_CHAR
                      IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_CHAR(J)))THEN
                         VALID = .TRUE.
                      ENDIF
                   ENDDO
                   DO J=1,NKEY_INT
                      IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_INT(J)))THEN
                         VALID = .TRUE.
                      ENDIF
                   ENDDO
                   DO J=1,NKEY_RE
                      IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_RE(J)))THEN
                         VALID = .TRUE.
                      ENDIF
                   ENDDO
                   DO J=1,NKEY_LOG
                      IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_LOG(J)))THEN
                         VALID = .TRUE.
                      ENDIF
                   ENDDO
                   IF(TRIM(VECT(I,K)).EQ.TRIM(STARTSTOP(2)))THEN
                      STOPL=.TRUE.
                   ENDIF
                   IF(.NOT.VALID.AND..NOT.STOPL)CALL CHECK_VALID(VECT(I,K))
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    STOPL = .FALSE.
    DO I=READMINI,READMAXI  !We search for the character keys
       IF(STOPL)EXIT
       DO K=1,READMAXJ
          IF(STOPL)EXIT
          IF(VECT(I,K).NE.'                    ')THEN
             IF(START)THEN !If we have a start key:
                IF(READMAXJ*(I-1)+K .GE.READMAXJ*(READMINI-1)+STARTATJ) THEN !If the position is beyond the start key:
                   IF(TRIM(VECT(I,K)).NE.'}')THEN  !If we don't have a stop key:
                      DO J=1,NKEY_CHAR
                         IF(ADJUSTL(TRIM(VECT(I,K))).EQ.ADJUSTL(TRIM(KEYVECTOR_CHAR(J))))THEN
                            VALVECTOR_CHAR(J)=ADJUSTL(TRIM(VECT(I,K+1)))
                            CHECKMISSING_CHAR(J) = .TRUE.
                         ENDIF
                      END DO
                   ELSE
                      STOPL = .TRUE.
                   ENDIF
                ENDIF
             ELSE  !If we don't have a start key:
                DO J=1,NKEY_CHAR
                   IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_CHAR(J)))THEN
                      VALVECTOR_CHAR(J)=TRIM(VECT(I,K+1))
                      CHECKMISSING_CHAR(J) = .TRUE.
                   ENDIF
                END DO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    STOPL = .FALSE.
    DO I=READMINI,READMAXI  !We search for the integer keys
       IF(STOPL)EXIT
       DO K=1,READMAXJ
          IF(STOPL)EXIT
          IF(VECT(I,K).NE.'                    ')THEN
             IF(START)THEN
                IF(READMAXJ*(I-1)+K .GE.READMAXJ*(READMINI-1)+STARTATJ) THEN
                   IF(ADJUSTL(TRIM(VECT(I,K))).NE.'}')THEN
                      DO J=1,NKEY_INT
                         IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_INT(J)))THEN
                            READ(VECT(I,K+1),*)VALVECTOR_INT(J)
                            CHECKMISSING_INT(J) = .TRUE.
                         ENDIF
                      ENDDO
                   ELSE
                      STOPL = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO J=1,NKEY_INT
                   IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_INT(J)))THEN
                      READ(VECT(I,K+1),*)VALVECTOR_INT(J)
                      CHECKMISSING_INT(J) = .TRUE.
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    STOPL = .FALSE.
    DO I=READMINI,READMAXI  !We search for the real keys
       IF(STOPL)EXIT
       DO K=1,READMAXJ
          IF(STOPL)EXIT
          IF(VECT(I,K).NE.'                    ')THEN
             IF(START)THEN
                IF(READMAXJ*(I-1)+K .GE.READMAXJ*(READMINI-1)+STARTATJ) THEN
                   IF(TRIM(VECT(I,K)).NE.'}')THEN
                      DO J=1,NKEY_RE
                         IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_RE(J)))THEN
                            READ(VECT(I,K+1),*)VALVECTOR_RE(J)
                            CHECKMISSING_RE(J) = .TRUE.
                         ENDIF
                      ENDDO
                   ELSE
                      STOPL = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO J=1,NKEY_RE
                   IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_RE(J)))THEN
                      READ(VECT(I,K+1),*)VALVECTOR_RE(J)
                      CHECKMISSING_RE(J) = .TRUE.
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    STOPL = .FALSE.
    DO I=1,READMAXI  !We search for the logical keys
       IF(STOPL)EXIT
       DO K=1,READMAXJ
          IF(STOPL)EXIT
          IF(VECT(I,K).NE.'                    ')THEN
             IF(START)THEN
                IF(READMAXJ*(I-1)+K .GE.READMAXJ*(READMINI-1)+STARTATJ) THEN
                   IF(TRIM(VECT(I,K)).NE.'}')THEN
                      DO J=1,NKEY_LOG
                         IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_LOG(J)))THEN
                            READ(VECT(I,K+1),*)VALVECTOR_LOG(J)
                            CHECKMISSING_LOG(J) = .TRUE.
                         END IF
                      END DO
                   ELSE
                      STOPL = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO J=1,NKEY_LOG
                   IF(TRIM(VECT(I,K)).EQ.TRIM(KEYVECTOR_LOG(J)))THEN
                      READ(VECT(I,K+1),*)VALVECTOR_LOG(J)
                      CHECKMISSING_LOG(J) = .TRUE.
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    !Check for missing keywords
    WRITE(*,*)' '
    DO I = 1,NKEY_CHAR
       IF(DEFAULTNONE .EQV..TRUE.)THEN
          IF(CHECKMISSING_CHAR(I).NEQV..TRUE..AND.TRIM(KEYVECTOR_CHAR(I)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(KEYVECTOR_CHAR(I)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',VALVECTOR_CHAR(I)
             STOP
          ENDIF
       ENDIF
       IF((CHECKMISSING_CHAR(I).NEQV..TRUE.))  WRITE(*,*)'# WARNING: variable ',TRIM(KEYVECTOR_CHAR(I)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO I = 1,NKEY_INT
       IF(DEFAULTNONE .EQV..TRUE.)THEN
          IF(CHECKMISSING_INT(I).NEQV..TRUE..AND.TRIM(KEYVECTOR_INT(I)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(KEYVECTOR_INT(I)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',VALVECTOR_INT(I)
             STOP
          ENDIF
       ENDIF
       IF((CHECKMISSING_INT(I).NEQV..TRUE.)) WRITE(*,*)'# WARNING: variable ',TRIM(KEYVECTOR_INT(I)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO I = 1,NKEY_RE
       IF(DEFAULTNONE .EQV..TRUE.)THEN
          IF(CHECKMISSING_RE(I).NEQV..TRUE..AND.TRIM(KEYVECTOR_RE(I)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(KEYVECTOR_RE(I)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',VALVECTOR_RE(I)
             STOP
          ENDIF
       ENDIF
       IF((CHECKMISSING_RE(I).NEQV..TRUE.)) WRITE(*,*)'# WARNING: variable ',TRIM(KEYVECTOR_RE(I)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO I = 1,NKEY_LOG
       IF(DEFAULTNONE .EQV..TRUE.)THEN
          IF(CHECKMISSING_LOG(I).NEQV..TRUE..AND.TRIM(KEYVECTOR_LOG(I)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(KEYVECTOR_LOG(I)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',VALVECTOR_LOG(I)
             STOP
          ENDIF
       ENDIF
       IF((CHECKMISSING_LOG(I).NEQV..TRUE.)) WRITE(*,*)'# WARNING: variable ',TRIM(KEYVECTOR_LOG(I)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    WRITE(*,*)' '

    DEALLOCATE(CHECKMISSING_CHAR,CHECKMISSING_RE, CHECKMISSING_INT, CHECKMISSING_LOG)

    ! Only rank 0 prints parameters if compiled with MPI
    WRITE(*,*)' '
    !    if (printRank() .eq. 1) then

    WRITE(*,*)"############### PARAMETERS USED FOR THIS RUN ################"
    IF(START)WRITE(*,*)"#  ",STARTSTOP(1)
    DO J=1,NKEY_INT
       WRITE(*,*)"#    ",TRIM(KEYVECTOR_INT(J)),VALVECTOR_INT(J)
    ENDDO

    DO J=1,NKEY_RE
       WRITE(*,*)"#    ",TRIM(KEYVECTOR_RE(J)),VALVECTOR_RE(J)
    ENDDO

    DO J=1,NKEY_CHAR
       WRITE(*,*)"#    ",TRIM(KEYVECTOR_CHAR(J)),VALVECTOR_CHAR(J)
    ENDDO

    DO J=1,NKEY_LOG
       WRITE(*,*)"#    ",TRIM(KEYVECTOR_LOG(J)),VALVECTOR_LOG(J)
    ENDDO
    IF(START)WRITE(*,*)"#  ",STARTSTOP(2)

    WRITE(*,*)' '

    !   endif

    IF(STOPPARSING)THEN
       WRITE(*,*)"" ; WRITE(*,*)"STOP key found. Stop parsing ... "; WRITE(*,*)""
       STOP
    ENDIF

    DEALLOCATE(VECT)

  END SUBROUTINE PARSING_KERNEL

  !> Check for valid keywords (checks for an = sign)
  !! \param invalidc Keyword to check.
  !!
  SUBROUTINE CHECK_VALID(INVALIDC)
    IMPLICIT NONE
    CHARACTER(1), ALLOCATABLE     ::  TEMPC(:)
    CHARACTER(LEN=*), INTENT(IN)  ::  INVALIDC
    CHARACTER(LEN=100)            ::  TEMPCFLEX
    INTEGER                       ::  L, LENC

    LENC=LEN(ADJUSTL(TRIM(INVALIDC)))
    IF(.NOT.ALLOCATED(TEMPC))ALLOCATE(TEMPC(LENC))
    DO L = 1,LEN(ADJUSTL(TRIM(INVALIDC)))
       TEMPCFLEX = ADJUSTL(TRIM(INVALIDC))
       TEMPC(L) = TEMPCFLEX(L:L)
       IF(TEMPC(L).EQ."=".AND.TEMPC(1).NE."#")THEN
          WRITE(*,*)"Input ERROR: ",ADJUSTL(TRIM(INVALIDC))," is not a valid keyword"
          STOP
       ENDIF
    ENDDO

  END SUBROUTINE CHECK_VALID

END MODULE KERNELPARSER_MOD
