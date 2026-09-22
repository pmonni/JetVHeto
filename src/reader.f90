!!======================================================================
!!! module to help read parameters and tabulated data
!!! (e.g. fixed-order distributions) from user-supplied card files.
module reader
  use types; use consts_dp
  use warnings_and_errors
  use sub_defs_io
  implicit none

  private
  public :: read_file_into_array
  public :: reader_opts
  public :: dble_val, int_val, string_val, is_present
  public :: assert_all_opts_used

  integer, parameter :: reader_nopts_max = 300, reader_optname_maxlen = 20, reader_optvalue_maxlen = 400
  type reader_opts
     character(len=reader_optvalue_maxlen) :: optvalues(reader_nopts_max), filename
     character(len=reader_optname_maxlen)  :: optnames(reader_nopts_max)
     logical :: used(reader_nopts_max)
     integer :: nopts
  end type reader_opts
  

contains

  !======================================================================
  ! This subroutines processes a card file: lines of the form
  !
  !   variable = value   # some comment
  !
  ! are read in and stored in the "opts" object, for later access
  ! through the dble_val(...), int_val(...) and string_val(...)
  ! functions. Variable names are treated case insensitive.
  !
  ! Lines that start with numeric information are assumed to be
  ! tabulated data. These are read into the used-supplied "array"
  ! variable, which on return is of exactly the right size to contain
  ! all the tabulated data.
  !
  ! Lines starting with the "#" or "!" characters are taken to be
  ! comments and ignored.
  !
  ! If the "max_col1" argument is supplied, reading stops at the first
  ! occurrence of tabulated data whose first column is greater than
  ! max_col1.
  subroutine read_file_into_array(filename, array, opts, max_col1)
    character(len=*),  intent(in)                 :: filename
    real(dp),          intent(inout), allocatable :: array(:,:)
    type(reader_opts), intent(out)                :: opts
    real(dp),          intent(in), optional       :: max_col1
    !---------------------------------------
    integer             :: ncols, nrows, iunit
    character(len=1000) :: line, copy
    integer :: ios, i
    logical :: first_row_done
    
    iunit = get_new_device()
    opts%filename = filename
    open (iunit, file=filename, status="OLD")
    opts%nopts = 0
    opts%used = .false.
    first_row_done = .false.
    nrows = 0
    do
      read(iunit,'(a)',iostat=ios) line
      if (ios /= 0) exit

      line = adjustl(line)
      
      if (line(1:1) == "#" .or. line(1:1) == "!") then
         ! we ignore comment lines
         cycle
      else if ((lge(line(1:1),'a') .and. lle(line(1:1),'z')) .or.&
           &   (lge(line(1:1),'A') .and. lle(line(1:1),'Z'))) then
         ! lines starting with a character indicate an option
         i = index(line,'=')
         if (i == 0) call wae_error('could not interpret the option line: '//trim(line))
         opts%nopts = opts%nopts + 1
         read(line(:i-1),*) opts%optnames(opts%nopts)
         ! NB: read with "*" formating takes in only the first word of 
         !     a string (any subsequent words, and any comments, are not read)
         read(line(i+1:),*) opts%optvalues(opts%nopts)
         call lowercase(opts%optnames(opts%nopts))
         cycle
      end if
      if (len_trim(line) == 0) cycle


      if (nrows == 0) then
        ! try to count the number of columns (NB: not very efficient...)
        ncols = 1
        copy = line
        do
          i = index(trim(copy),' ')
          if (i == 0) exit
          copy = adjustl(copy(i:))
          ncols = ncols + 1
        end do
        allocate(array(ncols,10)) ! start off with 10 rows
      end if
      
      if (nrows >= size(array,2)) call resize(array, nrows, 2*nrows)
      nrows = nrows + 1
      !write(0,*) trim(line)
      read(line,*) array(:,nrows)
      if (array(1,nrows) > max_col1) then
         nrows = nrows - 1
         exit
      end if
    end do
    
    call resize(array, nrows, nrows)
    close (iunit)
  end subroutine read_file_into_array

 
  !======================================================================
  !! Given the "opts" object, initialised by
  !! read_file_into_array(...), return the value of the option with
  !! name opt_name, interpreted as a double; if it is not present and
  !! the "default" argument was supplied, return that instead.
  !!
  !! (If it is not present and the default argument was not supplied,
  !! it gives an error).
  function dble_val(opts, opt_name, default) result(res)
    type(reader_opts), intent(inout)  :: opts
    character(len=*),  intent(in)  :: opt_name
    real(dp), intent(in), optional ::  default
    real(dp) :: res
    !----------------------
    integer :: i
    i = iopt(opts, opt_name)
    if (i == 0) then
       if (present(default)) then
          res = default
       else
          call wae_error('option processing','Could not find option '//opt_name)
       end if
    else
       read(opts%optvalues(i), *) res
    end if
  end function dble_val

  !======================================================================
  !! like dble_val, but interprets the option as an integer value
  function int_val(opts, opt_name, default) result(res)
    type(reader_opts), intent(inout)  :: opts
    character(len=*),  intent(in)  :: opt_name
    integer, intent(in), optional ::  default
    integer :: res
    !----------------------
    integer :: i
    i = iopt(opts, opt_name)
    if (i == 0) then
       if (present(default)) then
          res = default
       else
          call wae_error('option processing','Could not find option '//opt_name)
       end if
    else
       read(opts%optvalues(i), *) res
    end if
  end function int_val

  !======================================================================
  !! like dble_val, but interprets the option as a string: not only
  !! the first word of the string is taken, and the string should not
  !! be in quotes. (If you're trying to read a filename with spaces in
  !! it, you're in trouble!)
  function string_val(opts, opt_name, default) result(res)
    type(reader_opts), intent(inout)  :: opts
    character(len=*),  intent(in)  :: opt_name
    character(len=*), intent(in), optional ::  default
    character(len=reader_optvalue_maxlen) :: res
    !----------------------
    integer :: i
    i = iopt(opts, opt_name)
    if (i == 0) then
       if (present(default)) then
          res = default
       else
          call wae_error('option processing','Could not find option '//opt_name)
       end if
    else
       read(opts%optvalues(i), *) res
    end if
  end function string_val
  
  !======================================================================
  !! return "true" if the option is present
  logical function is_present(opts, opt_name) result(res)
    type(reader_opts), intent(inout)  :: opts
    character(len=*),  intent(in)  :: opt_name
    integer :: i
    i = iopt(opts, opt_name)
    res = .not. (i == 0)
  end function is_present

  !======================================================================
  !! stops the program if there were options in the user's card that
  !! had not been read in.
  subroutine assert_all_opts_used(opts)
    type(reader_opts), intent(in) :: opts
    integer i, unused
    unused = 0
    do i = 1, opts%nopts
       if (.not. opts%used(i)) then
          unused = unused + 1
          write(0,*) 'Error: option '//trim(opts%optnames(i))//&
               &     " in file "//trim(opts%filename)//" was not recognized"
       end if
    end do
    if (unused /= 0) stop
  end subroutine assert_all_opts_used
  

  !======================================================================
  !! for internal use only
  integer function iopt(opts, opt_name)
    type(reader_opts), intent(inout) :: opts
    character(len=*),  intent(in) :: opt_name
    !----------------------------------------
    character(len=len(opt_name)) :: lc_opt_name
    lc_opt_name = opt_name
    call lowercase(lc_opt_name)
    do iopt = 1, opts%nopts
       if (trim(lc_opt_name) == trim(opts%optnames(iopt))) exit
    end do
    if (iopt > opts%nopts) then
       iopt = 0
    else
       opts%used(iopt) = .true.
    end if
  end function iopt
  


  !======================================================================
  !! resize an array so that the numbers of number of rows is changed
  !! to nrows_wanted
  !
  !! intended for internal use only
  subroutine resize(array, nrows_used, nrows_wanted)
    real(dp), intent(inout), allocatable :: array(:,:)
    integer,  intent(in) :: nrows_used, nrows_wanted
    !------------------------------
    real(dp), allocatable :: array_tmp(:,:)
    integer  :: ncols, nrows_min

    nrows_min = min(nrows_used, nrows_wanted)
    ncols = size(array,1)

    allocate(array_tmp(ncols, nrows_min))
    array_tmp = array(:,:nrows_min)

    deallocate(array)
    allocate(array(ncols, nrows_wanted))
    array(:,:nrows_min) = array_tmp(:,:nrows_min)

    deallocate(array_tmp)
  end subroutine resize
  

  !----------------------------------------------------------------------
  !! routine inspired by from http://www.davidgsimpson.com/software/chcase_f90.txt to change case
  subroutine lowercase(string)
    implicit none
    character(len=*), intent(inout) :: string
    integer :: i, delta
    delta = iachar('a') - iachar('A')
    do i = 1, len_trim(string)
       if (lge(string(i:i),'A') .and. lle(string(i:i),'Z')) then
          string(i:i) = achar(iachar(string(i:i)) + delta)
       end if
    end do
  end subroutine lowercase

 
end module reader

! ! temporary test program
! program test_reader
!   use reader
!   use types
!   implicit none
!   type(reader_opts) :: opts
!   real(dp), allocatable :: array(:,:)
!   character(len=30) :: M
! 
!   call read_file_into_array("/dev/stdin",array,opts)
!   M = string_val(opts,"m","30")
!   write(0,*) "M = ", M
!   call assert_all_opts_used(opts)
! end program test_reader
