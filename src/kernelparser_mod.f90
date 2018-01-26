!> Some general parsing functions.
!! \ingroup PROGRESS
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
MODULE kernelparser_mod

  USE openfiles_mod
  !  use parallel_mod

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: dp = KIND(1.0d0)

  PUBLIC :: parsing_kernel

CONTAINS

  !> The general parsing function.
  !! It is used to vectorize a set of "keywords" "value" pairs
  !! as included in a general input file.
  !! \note This parsing strategy can only parse a file of
  !! 500 lines by 500 words.
  !! \warning If the length of variable vect is changed, this could produce a
  !! segmentation fault.
  !!
  SUBROUTINE parsing_kernel(keyvector_char,valvector_char&
       ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
       keyvector_log,valvector_log,filename,startstop)
    IMPLICIT NONE
    CHARACTER(1), ALLOCATABLE              ::  tempc(:)
    CHARACTER(100), ALLOCATABLE            ::  vect(:,:)
    CHARACTER(50)                          ::  keyvector_char(:), keyvector_int(:), keyvector_log(:), keyvector_re(:)
    CHARACTER(100)                         ::  valvector_char(:)
    CHARACTER(len=*)                       ::  filename
    CHARACTER(len=*), INTENT(in), OPTIONAL ::  startstop(2)
    CHARACTER(len=100)                      ::  tempcflex
    INTEGER                                ::  i, io_control, ios, j
    INTEGER                                ::  k, l, lenc, nkey_char
    INTEGER                                ::  nkey_int, nkey_log, nkey_re, readmaxi
    INTEGER                                ::  readmaxj, readmini, valvector_int(:)
    INTEGER                                ::  startatj, totalwords
    LOGICAL                                ::  start, stopl, valid, valvector_log(:), stopparsing, defaultnone
    LOGICAL, ALLOCATABLE                   ::  checkmissing_char(:), checkmissing_int(:), checkmissing_log(:), checkmissing_re(:)
    REAL(dp)                               ::  valvector_re(:)

    readmaxi = 5 ; readmaxj = 1000
    ALLOCATE(vect(readmaxi,readmaxj))
    nkey_char = SIZE(keyvector_char,dim=1)
    nkey_re = SIZE(keyvector_re,dim=1)
    nkey_int = SIZE(keyvector_int,dim=1)
    nkey_log = SIZE(keyvector_log,dim=1)

    CALL open_file_to_read(io_control,filename)

    ALLOCATE(checkmissing_char(nkey_char),checkmissing_re(nkey_re), &
         checkmissing_int(nkey_int), checkmissing_log(nkey_log))

    !Initialize the checkmissing flags and the vect arrays
    checkmissing_char = .FALSE.
    checkmissing_re = .FALSE.
    checkmissing_int = .FALSE.
    checkmissing_log = .FALSE.
    stopparsing = .FALSE.
    defaultnone = .FALSE.
    vect = '                    '

    DO i=1,readmaxi !Here we read all the input into vect
       READ(io_control,*,iostat=ios)(vect(i,j),j=1,readmaxj)
    END DO

    CLOSE(io_control)

    !Look up for floating hashes (#)
    totalwords = 0
    DO i=1,readmaxi
       DO k=1,readmaxj
          IF(ADJUSTL(TRIM(vect(i,k))).NE."")totalwords = totalwords + 1
          IF(ADJUSTL(TRIM(vect(i,k))).EQ."#")THEN
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
          IF(ADJUSTL(TRIM(vect(i,k))).EQ."STOP{}")stopparsing = .TRUE.
          IF(ADJUSTL(TRIM(vect(i,k))).EQ."DEFAULTNONE")defaultnone = .TRUE.
       ENDDO
    ENDDO

    IF(totalwords > readmaxi*readmaxj - 100) THEN
       WRITE(*,*)""; WRITE(*,*)"Stopping ... Maximum allowed (keys + values + comments) words close to the limit "
       WRITE(*,*)"Increase the readmaxj variable in the parsing_kernel subroutine or reduce the comments in the input"
       STOP
    ENDIF

    !Look up for boundaries
    readmini=1
    start=.FALSE.
    IF(PRESENT(startstop))THEN
       DO i=1,readmaxi
          DO k=1,readmaxj
             IF(TRIM(vect(i,k)).EQ.TRIM(startstop(1)))THEN
                readmini=i
                startatj=k
                start=.TRUE.
             ENDIF
             IF(start.AND.TRIM(vect(i,k)).EQ.TRIM(startstop(2)))THEN
                readmaxi=i
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    WRITE(*,*)startstop

    ! Look for invalid characters if startstop is present
    IF(start)THEN
       start=.FALSE.
       stopl=.FALSE.
       DO i=readmini,readmaxi
          DO k=1,readmaxj
             IF(TRIM(vect(i,k)).EQ.TRIM(startstop(1)))start=.TRUE.
             valid = .FALSE.
             IF(start)THEN
                IF(vect(i,k).NE.'                    ')THEN
                   DO j=1,nkey_char
                      IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_char(j)))THEN
                         valid = .TRUE.
                      ENDIF
                   ENDDO
                   DO j=1,nkey_int
                      IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_int(j)))THEN
                         valid = .TRUE.
                      ENDIF
                   ENDDO
                   DO j=1,nkey_re
                      IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_re(j)))THEN
                         valid = .TRUE.
                      ENDIF
                   ENDDO
                   DO j=1,nkey_log
                      IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_log(j)))THEN
                         valid = .TRUE.
                      ENDIF
                   ENDDO
                   IF(TRIM(vect(i,k)).EQ.TRIM(startstop(2)))THEN
                      stopl=.TRUE.
                   ENDIF
                   IF(.NOT.valid.AND..NOT.stopl)CALL check_valid(vect(i,k))
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    stopl = .FALSE.
    DO i=readmini,readmaxi  !We search for the character keys
       IF(stopl)EXIT
       DO k=1,readmaxj
          IF(stopl)EXIT
          IF(vect(i,k).NE.'                    ')THEN
             IF(start)THEN !If we have a start key:
                IF(readmaxj*(i-1)+k .GE.readmaxj*(readmini-1)+startatj) THEN !If the position is beyond the start key:
                   IF(TRIM(vect(i,k)).NE.'}')THEN  !If we don't have a stop key:
                      DO j=1,nkey_char
                         IF(ADJUSTL(TRIM(vect(i,k))).EQ.ADJUSTL(TRIM(keyvector_char(j))))THEN
                            valvector_char(j)=ADJUSTL(TRIM(vect(i,k+1)))
                            checkmissing_char(j) = .TRUE.
                         ENDIF
                      END DO
                   ELSE
                      stopl = .TRUE.
                   ENDIF
                ENDIF
             ELSE  !If we don't have a start key:
                DO j=1,nkey_char
                   IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_char(j)))THEN
                      valvector_char(j)=TRIM(vect(i,k+1))
                      checkmissing_char(j) = .TRUE.
                   ENDIF
                END DO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    stopl = .FALSE.
    DO i=readmini,readmaxi  !We search for the integer keys
       IF(stopl)EXIT
       DO k=1,readmaxj
          IF(stopl)EXIT
          IF(vect(i,k).NE.'                    ')THEN
             IF(start)THEN
                IF(readmaxj*(i-1)+k .GE.readmaxj*(readmini-1)+startatj) THEN
                   IF(ADJUSTL(TRIM(vect(i,k))).NE.'}')THEN
                      DO j=1,nkey_int
                         IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_int(j)))THEN
                            READ(vect(i,k+1),*)valvector_int(j)
                            checkmissing_int(j) = .TRUE.
                         ENDIF
                      ENDDO
                   ELSE
                      stopl = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO j=1,nkey_int
                   IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_int(j)))THEN
                      READ(vect(i,k+1),*)valvector_int(j)
                      checkmissing_int(j) = .TRUE.
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    stopl = .FALSE.
    DO i=readmini,readmaxi  !We search for the real keys
       IF(stopl)EXIT
       DO k=1,readmaxj
          IF(stopl)EXIT
          IF(vect(i,k).NE.'                    ')THEN
             IF(start)THEN
                IF(readmaxj*(i-1)+k .GE.readmaxj*(readmini-1)+startatj) THEN
                   IF(TRIM(vect(i,k)).NE.'}')THEN
                      DO j=1,nkey_re
                         IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_re(j)))THEN
                            READ(vect(i,k+1),*)valvector_re(j)
                            checkmissing_re(j) = .TRUE.
                         ENDIF
                      ENDDO
                   ELSE
                      stopl = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO j=1,nkey_re
                   IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_re(j)))THEN
                      READ(vect(i,k+1),*)valvector_re(j)
                      checkmissing_re(j) = .TRUE.
                   ENDIF
                ENDDO
             ENDIF
          ELSE
             EXIT
          ENDIF
       ENDDO
    ENDDO

    stopl = .FALSE.
    DO i=1,readmaxi  !We search for the logical keys
       IF(stopl)EXIT
       DO k=1,readmaxj
          IF(stopl)EXIT
          IF(vect(i,k).NE.'                    ')THEN
             IF(start)THEN
                IF(readmaxj*(i-1)+k .GE.readmaxj*(readmini-1)+startatj) THEN
                   IF(TRIM(vect(i,k)).NE.'}')THEN
                      DO j=1,nkey_log
                         IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_log(j)))THEN
                            READ(vect(i,k+1),*)valvector_log(j)
                            checkmissing_log(j) = .TRUE.
                         END IF
                      END DO
                   ELSE
                      stopl = .TRUE.
                   ENDIF
                ENDIF
             ELSE
                DO j=1,nkey_log
                   IF(TRIM(vect(i,k)).EQ.TRIM(keyvector_log(j)))THEN
                      READ(vect(i,k+1),*)valvector_log(j)
                      checkmissing_log(j) = .TRUE.
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
    DO i = 1,nkey_char
       IF(defaultnone .EQV..TRUE.)THEN
          IF(checkmissing_char(i).NEQV..TRUE..AND.TRIM(keyvector_char(i)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(keyvector_char(i)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',valvector_char(i)
             STOP
          ENDIF
       ENDIF
       IF(checkmissing_char(i).NEQV..TRUE.)  WRITE(*,*)'WARNING: variable ',TRIM(keyvector_char(i)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO i = 1,nkey_int
       IF(defaultnone .EQV..TRUE.)THEN
          IF(checkmissing_int(i).NEQV..TRUE..AND.TRIM(keyvector_int(i)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(keyvector_int(i)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',valvector_int(i)
             STOP
          ENDIF
       ENDIF
       IF(checkmissing_int(i).NEQV..TRUE.) WRITE(*,*)'WARNING: variable ',TRIM(keyvector_int(i)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO i = 1,nkey_re
       IF(defaultnone .EQV..TRUE.)THEN
          IF(checkmissing_re(i).NEQV..TRUE..AND.TRIM(keyvector_re(i)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(keyvector_re(i)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',valvector_re(i)
             STOP
          ENDIF
       ENDIF
       IF(checkmissing_re(i).NEQV..TRUE.) WRITE(*,*)'WARNING: variable ',TRIM(keyvector_re(i)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    DO i = 1,nkey_log
       IF(defaultnone .EQV..TRUE.)THEN
          IF(checkmissing_log(i).NEQV..TRUE..AND.TRIM(keyvector_log(i)).NE."DUMMY=")THEN
             WRITE(*,*)'ERROR: variable ',TRIM(keyvector_log(i)),&
                  ' is missing. Set this variable or remove the DEFAULTNONE keyword from the input file...'
             WRITE(*,*)'Default value is:',valvector_log(i)
             STOP
          ENDIF
       ENDIF
       IF(checkmissing_log(i).NEQV..TRUE.) WRITE(*,*)'WARNING: variable ',TRIM(keyvector_log(i)),&
            ' is missing. I will use a default value instead ...'
    ENDDO
    WRITE(*,*)' '

    DEALLOCATE(checkmissing_char,checkmissing_re, checkmissing_int, checkmissing_log)

    ! Only rank 0 prints parameters if compiled with MPI
    WRITE(*,*)' '
    !    if (printRank() .eq. 1) then

    WRITE(*,*)"############### Parameters used for this run ################"
    IF(start)WRITE(*,*)" ",startstop(1)
    DO j=1,nkey_int
       WRITE(*,*)" ",TRIM(keyvector_int(j)),valvector_int(j)
    ENDDO

    DO j=1,nkey_re
       WRITE(*,*)" ",TRIM(keyvector_re(j)),valvector_re(j)
    ENDDO

    DO j=1,nkey_char
       WRITE(*,*)" ",TRIM(keyvector_char(j)),valvector_char(j)
    ENDDO

    DO j=1,nkey_log
       WRITE(*,*)" ",TRIM(keyvector_log(j)),valvector_log(j)
    ENDDO
    IF(start)WRITE(*,*)" ",startstop(2)

    !   endif
    WRITE(*,*)' '

    IF(stopparsing)THEN
       WRITE(*,*)"" ; WRITE(*,*)"STOP key found. Stop parsing ... "; WRITE(*,*)""
       STOP
    ENDIF

    DEALLOCATE(vect)

  END SUBROUTINE parsing_kernel

  !> Check for valid keywords (checks for an = sign)
  !! \param invalidc Keyword to check.
  !!
  SUBROUTINE check_valid(invalidc)
    IMPLICIT NONE
    CHARACTER(1), ALLOCATABLE     ::  tempc(:)
    CHARACTER(len=*), INTENT(in)  ::  invalidc
    CHARACTER(len=100)            ::  tempcflex
    INTEGER                       ::  l, lenc

    lenc=LEN(ADJUSTL(TRIM(invalidc)))
    IF(.NOT.ALLOCATED(tempc))ALLOCATE(tempc(lenc))
    DO l = 1,LEN(ADJUSTL(TRIM(invalidc)))
       tempcflex = ADJUSTL(TRIM(invalidc))
       tempc(l) = tempcflex(l:l)
       IF(tempc(l).EQ."=".AND.tempc(1).NE."#")THEN
          WRITE(*,*)"Input ERROR: ",ADJUSTL(TRIM(invalidc))," is not a valid keyword"
          STOP
       ENDIF
    ENDDO

  END SUBROUTINE check_valid


END MODULE kernelparser_mod
