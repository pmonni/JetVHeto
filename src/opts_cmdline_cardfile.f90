module opts_cmdline_cardfile
  use sub_defs_io
  use reader
  use types
  implicit none
  
  private 
  public :: dble_val_opts_or_cmdline, int_val_opts_or_cmdline, string_val_opts_or_cmdline
  public :: var_or_xVar 

contains
  ! reads an option named opt_name from the "opts" objects (which
  ! comes from a cardfile), or from the command line (it should then
  ! be -opt_name).
  !
  ! if neither is present, and default is supplied, return default,
  ! otherwise give an error.
  !
  ! If the override_source string is present, then if the option is
  ! set from opts or from cmdline, i.e. overriding the "default"
  ! value, then inform the user telling them that the value from
  ! override_source was overwritten.
  !
  ! If give_error is present and true, then the same option may not be
  ! supplied in both the "opts" object and the command line.
  function dble_val_opts_or_cmdline(opts, opt_name, default, give_error) result(res)
    type(reader_opts),   intent(inout) :: opts
    character(len=*),    intent(in) :: opt_name
    real(dp), optional,  intent(in) :: default
    logical, intent(in), optional :: give_error
    real(dp) :: res

    if (is_present(opts, opt_name)) then
       res = dble_val(opts, opt_name)
       ! see if the command-line overrides it
       if (log_val_opt("-"//opt_name)) then
          write(0,*) "**** Command line overrides file's value for ",trim(opt_name)
          if (present(give_error)) then
             if (give_error) then 
                write(0,*) "Overriding option "//trim(opt_name)//" is forbidden"
                stop 
             end if
          end if
          res = dble_val_opt("-"//opt_name)
       end if
    else
       res = dble_val_opt("-"//opt_name, default)
    end if

  end function dble_val_opts_or_cmdline

  function int_val_opts_or_cmdline(opts, opt_name, default, give_error) result(res)
    type(reader_opts),   intent(inout) :: opts
    character(len=*),    intent(in) :: opt_name
    integer, optional,  intent(in) :: default
    logical, intent(in), optional :: give_error
    integer :: res

    if (is_present(opts, opt_name)) then
       res = int_val(opts, opt_name)
       ! see if the command-line overrides it
       if (log_val_opt("-"//opt_name)) then
          write(0,*) "**** Command line overrides file's value for ",trim(opt_name)
          if (present(give_error)) then
             if (give_error) then 
                write(0,*) "Overriding option "//trim(opt_name)//" is forbidden"
                stop 
             end if
          end if
          res = int_val_opt("-"//opt_name)
       end if
    else
       res = int_val_opt("-"//opt_name, default)
    end if

  end function int_val_opts_or_cmdline


  function string_val_opts_or_cmdline(opts, opt_name, default, give_error) result(res)
    type(reader_opts),   intent(inout) :: opts
    character(len=*),    intent(in) :: opt_name
    character(len=*), optional,  intent(in) :: default
    logical, intent(in), optional :: give_error
    character(len=100) :: res

    if (is_present(opts, opt_name)) then
       res = string_val(opts, opt_name)
       ! see if the command-line overrides it
       if (log_val_opt("-"//opt_name)) then
          write(0,*) "**** Command line overrides file's value for ",trim(opt_name)
          if (present(give_error)) then
             if (give_error) then
                write(0,*) "Overriding option "//trim(opt_name)//" is forbidden"
                stop 
             end if
          end if
          res = string_val_opt("-"//opt_name)
       end if
    else
       res = string_val_opt("-"//opt_name, default)
    end if

  end function string_val_opts_or_cmdline


  !----------------------------------------------------------------------
  ! function to determine the value of a scale named varname. It can
  ! be provided in the fixed-order card-file or on the command line,
  ! both as absolute values and as a fraction of the boson mass (by
  ! preceding the name with "x").
  !
  ! override_allowed indicates if the command-line can override a value 
  ! already specified in the card file.
  function var_or_xVar(varname, xdefault, opts, M, override_allowed) result(var)
    character(len=*), intent(in)  :: varname
    real(dp),         intent(in)  :: xdefault,M
    type(reader_opts), intent(inout) :: opts 
    logical,          intent(in)  :: override_allowed
    real(dp)                      :: var
    !--------------------------------------
    real(dp) :: xvar

    if (log_val_opt("-"//varname) .and. log_val_opt("-x"//varname)) then
       write(0,*) "Error: options ","-"//varname, " and ", "-x"//varname
       write(0,*) "were both present on the command line; this is not allowed"
       stop
    end if
    
    if (is_present(opts,varname) .and. is_present(opts,"x"//varname)) then
       write(0,*) "Error: options ",varname, " and ", "x"//varname
       write(0,*) "were both present in the input file; this is not allowed"
       stop
    end if

    ! try getting a value from the card file
    xvar = dble_val(opts,"x"//varname,xdefault)
    var  = dble_val(opts,varname,xvar*M)

    ! then let the command-line override it
    xvar = dble_val_opt("-x"//varname,var/M)
    var  = dble_val_opt("-"//varname,xvar*M)

    ! give the user some warnings / errors if things are being overridden
    if (      (log_val_opt("-"//varname) .or. log_val_opt("-x"//varname))&
       &.and. (is_present(opts,varname)  .or. is_present(opts,"x"//varname))) then
       if (override_allowed) then
          write(0,*) "******************* WARNING: ", varname, " (or ", "x"//varname, ") were specified"
          write(0,*) "both in cardfile and command-line; command-line value of ",varname,"=",var
          write(0,*) "will be used"
       else
          write(0,*) "ERROR: card file and command-line both specified values for"
          write(0,*) varname, " (or ", "x"//varname, "). This is not allowed."
          stop
       end if
    end if
    
  end function var_or_xVar
  
    


end module opts_cmdline_cardfile
