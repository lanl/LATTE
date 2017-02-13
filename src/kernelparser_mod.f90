!> Some general parsing functions.
!! \ingroup PROGRESS
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
module kernelparser_mod

  use openfiles_mod
  use parallel_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: parsing_kernel

contains

  !> The general parsing function.
  !! It is used to vectorize a set of "keywords" "value" pairs
  !! as included in a general input file.
  !! \note This parsing strategy can only parse a file of
  !! 500 lines by 500 words.
  !! \warning If the length of variable vect is changed, this could produce a
  !! segmentation fault.
  !!
  subroutine parsing_kernel(keyvector_char,valvector_char&
      ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
      keyvector_log,valvector_log,filename,startstop)
    implicit none
    character(1), allocatable              ::  tempc(:)
    character(100), allocatable            ::  vect(:,:)
    character(50)                          ::  keyvector_char(:), keyvector_int(:), keyvector_log(:), keyvector_re(:)
    character(100)                         ::  valvector_char(:)
    character(len=*)                       ::  filename
    character(len=*), intent(in), optional ::  startstop(2)
    character(len=100)                      ::  tempcflex
    integer                                ::  i, io_control, ios, j
    integer                                ::  k, l, lenc, nkey_char
    integer                                ::  nkey_int, nkey_log, nkey_re, readmaxi
    integer                                ::  readmaxj, readmini, valvector_int(:)
    integer                                ::  startatj, totalwords
    logical                                ::  start, stopl, valid, valvector_log(:), stopparsing
    logical, allocatable                   ::  checkmissing_char(:), checkmissing_int(:), checkmissing_log(:), checkmissing_re(:)
    real(dp)                               ::  valvector_re(:)

    readmaxi = 5 ; readmaxj = 1000
    allocate(vect(readmaxi,readmaxj))
    nkey_char = size(keyvector_char,dim=1)
    nkey_re = size(keyvector_re,dim=1)
    nkey_int = size(keyvector_int,dim=1)
    nkey_log = size(keyvector_log,dim=1)

    call open_file_to_read(io_control,filename)

    allocate(checkmissing_char(nkey_char),checkmissing_re(nkey_re), &
      checkmissing_int(nkey_int), checkmissing_log(nkey_log))

    !Initialize the checkmissing flags and the vect array
    checkmissing_char = .false.
    checkmissing_re = .false.
    checkmissing_int = .false.
    checkmissing_log = .false.
    stopparsing = .false.
    vect = '                    '

    do i=1,readmaxi !Here we read all the input into vect
      read(io_control,*,iostat=ios)(vect(i,j),j=1,readmaxj)
    end do

    close(io_control)

    !Look up for floating dashes (#)
    totalwords = 0
    do i=1,readmaxi
      do k=1,readmaxj
        if(adjustl(trim(vect(i,k))).ne."")totalwords = totalwords + 1
        if(adjustl(trim(vect(i,k))).eq."#")then
          write(*,*)"ERROR in the the input file ..."
          write(*,*)" "
          write(*,*)"For this parsing routine everything is a comment by default unless theres an = sign"
          write(*,*)"next to a word in which case it will be recognized as a keyword."
          write(*,*)"This parser does not accept floating dashes (# symbol) to make sure"
          write(*,*)"that we are sure that we want to comment a specific keyword"
          write(*,*)" "
          write(*,*)"Make sure all the commented keywords have a # right next them"
          write(*,*)"   "
          write(*,*)"   The following commented keyword is correct: "
          write(*,*)"                #KeyWord= 1 "
          write(*,*)" "
          write(*,*)"   The following commented keyword is NOT correct: "
          write(*,*)"                # KeyWord= 1 "
          write(*,*)" "
          stop
        endif
        if(adjustl(trim(vect(i,k))).eq."STOP{}")stopparsing = .true.
      enddo
    enddo

     if(totalwords > readmaxi*readmaxj - 100) then
      write(*,*)""; write(*,*)"Stopping ... Maximum allowed (keys + values + comments) words close to the limit "
      write(*,*)"Increase the readmaxj variable in the parsing_kernel subroutine or reduce the comments in the input"
      stop
     endif

    !Look up for boundaries
    readmini=1
    start=.false.
    if(present(startstop))then
      do i=1,readmaxi
        do k=1,readmaxj
          if(trim(vect(i,k)).eq.trim(startstop(1)))then
            readmini=i
            startatj=k
            start=.true.
          endif
          if(start.and.trim(vect(i,k)).eq.trim(startstop(2)))then
            readmaxi=i
          endif
        enddo
      enddo
    endif
    write(*,*)start,startstop

    ! Look for invalid characters if startstop is present
    if(start)then
      start=.false.
      stopl=.false.
      do i=readmini,readmaxi
        do k=1,readmaxj
          if(trim(vect(i,k)).eq.trim(startstop(1)))start=.true.
          valid = .false.
          if(start)then
            if(vect(i,k).ne.'                    ')then
              do j=1,nkey_char
                if(trim(vect(i,k)).eq.trim(keyvector_char(j)))then
                  valid = .true.
                endif
              enddo
              do j=1,nkey_int
                if(trim(vect(i,k)).eq.trim(keyvector_int(j)))then
                  valid = .true.
                endif
              enddo
              do j=1,nkey_re
                if(trim(vect(i,k)).eq.trim(keyvector_re(j)))then
                  valid = .true.
                endif
              enddo
              do j=1,nkey_log
                if(trim(vect(i,k)).eq.trim(keyvector_log(j)))then
                  valid = .true.
                endif
              enddo
              if(trim(vect(i,k)).eq.trim(startstop(2)))then
                stopl=.true.
              endif
              if(.not.valid.and..not.stopl)call check_valid(vect(i,k))
            endif
          endif
        enddo
      enddo
    endif

    stopl = .false.
    do i=readmini,readmaxi  !We search for the character keys
      if(stopl)exit
      do k=1,readmaxj
        if(stopl)exit
        if(vect(i,k).ne.'                    ')then
          if(start)then !If we have a start key:
            if(readmaxj*(i-1)+k .ge.readmaxj*(readmini-1)+startatj) then !If the position is beyond the start key:
              if(trim(vect(i,k)).ne.'}')then  !If we don't have a stop key:
                do j=1,nkey_char
                  if(adjustl(trim(vect(i,k))).eq.adjustl(trim(keyvector_char(j))))then
                    valvector_char(j)=adjustl(trim(vect(i,k+1)))
                    checkmissing_char(j) = .true.
                  endif
                end do
              else
                stopl = .true.
              endif
            endif
          else  !If we don't have a start key:
            do j=1,nkey_char
              if(trim(vect(i,k)).eq.trim(keyvector_char(j)))then
                valvector_char(j)=trim(vect(i,k+1))
                 checkmissing_char(j) = .true.
              endif
            end do
          endif
        else
          exit
        end if
      end do
    end do

    stopl = .false.
    do i=readmini,readmaxi  !We search for the integer keys
      if(stopl)exit
      do k=1,readmaxj
        if(stopl)exit
        if(vect(i,k).ne.'                    ')then
          if(start)then
            if(readmaxj*(i-1)+k .ge.readmaxj*(readmini-1)+startatj) then
              if(adjustl(trim(vect(i,k))).ne.'}')then
                do j=1,nkey_int
                  if(trim(vect(i,k)).eq.trim(keyvector_int(j)))then
                    read(vect(i,k+1),*)valvector_int(j)
                    checkmissing_int(j) = .true.
                  end if
                end do
              else
                stopl = .true.
              endif
            endif
          else
            do j=1,nkey_int
              if(trim(vect(i,k)).eq.trim(keyvector_int(j)))then
                read(vect(i,k+1),*)valvector_int(j)
                checkmissing_int(j) = .true.
              end if
            end do
          endif
        else
          exit
        end if
      end do
    end do

    stopl = .false.
    do i=readmini,readmaxi  !We search for the real keys
      if(stopl)exit
      do k=1,readmaxj
        if(stopl)exit
        if(vect(i,k).ne.'                    ')then
          if(start)then
            if(readmaxj*(i-1)+k .ge.readmaxj*(readmini-1)+startatj) then
              if(trim(vect(i,k)).ne.'}')then
                do j=1,nkey_re
                  if(trim(vect(i,k)).eq.trim(keyvector_re(j)))then
                    read(vect(i,k+1),*)valvector_re(j)
                    checkmissing_re(j) = .true.
                  end if
                end do
              else
                stopl = .true.
              endif
            endif
          else
            do j=1,nkey_re
              if(trim(vect(i,k)).eq.trim(keyvector_re(j)))then
                read(vect(i,k+1),*)valvector_re(j)
                checkmissing_re(j) = .true.
              end if
            end do
          endif
        else
          exit
        end if
      end do
    end do

    stopl = .false.
    do i=1,readmaxi  !We search for the logical keys
      if(stopl)exit
      do k=1,readmaxj
        if(stopl)exit
        if(vect(i,k).ne.'                    ')then
          if(start)then
            if(readmaxj*(i-1)+k .ge.readmaxj*(readmini-1)+startatj) then
              if(trim(vect(i,k)).ne.'}')then
                do j=1,nkey_log
                  if(trim(vect(i,k)).eq.trim(keyvector_log(j)))then
                    read(vect(i,k+1),*)valvector_log(j)
                    checkmissing_log(j) = .true.
                  end if
                end do
              else
                stopl = .true.
              endif
            endif
          else
            do j=1,nkey_log
              if(trim(vect(i,k)).eq.trim(keyvector_log(j)))then
                read(vect(i,k+1),*)valvector_log(j)
                checkmissing_log(j) = .true.
              end if
            end do
          endif
        else
          exit
        end if
      end do
    end do

    !Check for missing keywords
    write(*,*)' '
    do i = 1,nkey_char
      if(checkmissing_char(i).neqv..true.) write(*,*)'WARNING: variable ',trim(keyvector_char(i)),&
        ' is missing. I will use a default value instead ...'
    end do
    do i = 1,nkey_int
      if(checkmissing_int(i).neqv..true.) write(*,*)'WARNING: variable ',trim(keyvector_int(i)),&
        ' is missing. I will use a default value instead ...'
    end do
    do i = 1,nkey_re
      if(checkmissing_re(i).neqv..true.) write(*,*)'WARNING: variable ',trim(keyvector_re(i)),&
        ' is missing. I will use a default value instead ...'
    end do
    do i = 1,nkey_log
      if(checkmissing_log(i).neqv..true.) write(*,*)'WARNING: variable ',trim(keyvector_log(i)),&
        ' is missing. I will use a default value instead ...'
    end do
    write(*,*)' '

    deallocate(checkmissing_char,checkmissing_re, checkmissing_int, checkmissing_log)

    ! Only rank 0 prints parameters
    write(*,*)' '
    if (printRank() .eq. 1) then

      write(*,*)"############### Parameters used for this run ################"
      if(start)write(*,*)" ",startstop(1)
      do j=1,nkey_int
        write(*,*)" ",trim(keyvector_int(j)),valvector_int(j)
      end do

      do j=1,nkey_re
        write(*,*)" ",trim(keyvector_re(j)),valvector_re(j)
      end do

      do j=1,nkey_char
        write(*,*)" ",trim(keyvector_char(j)),valvector_char(j)
      end do

      do j=1,nkey_log
        write(*,*)" ",trim(keyvector_log(j)),valvector_log(j)
      end do
      if(start)write(*,*)" ",startstop(2)

    endif
    write(*,*)' '

    if(stopparsing)then
      write(*,*)"" ; write(*,*)"STOP key found. Stop parsing ... "; write(*,*)""
      stop
    endif

    deallocate(vect)

  end subroutine parsing_kernel

  !> Check for valid keywords (checks for an = sign)
  !! \param invalidc Keyword to check.
  !!
  subroutine check_valid(invalidc)
    implicit none
    character(1), allocatable     ::  tempc(:)
    character(len=*), intent(in)  ::  invalidc
    character(len=100)            ::  tempcflex
    integer                       ::  l, lenc

    lenc=len(adjustl(trim(invalidc)))
    if(.not.allocated(tempc))allocate(tempc(lenc))
    do l = 1,len(adjustl(trim(invalidc)))
      tempcflex = adjustl(trim(invalidc))
      tempc(l) = tempcflex(l:l)
      if(tempc(l).eq."=".and.tempc(1).ne."#")then
        write(*,*)"Input ERROR: ",adjustl(trim(invalidc))," is not a valid keyword"
        stop
      endif
    enddo

  end subroutine check_valid


end module kernelparser_mod
