!! The coefficient functions for the resummation have been derived
!! from de Florian and Grazzini, hep-ph/0108273, Eqs.(46-47), with a
!! "hard" part that we put into H(1) separated out (this is so that
!! the hard part has scale alphas(muR), and the C1 part a scale
!! alphas(\sim ptveto)).
!!
!! The expressions used should correspond to Eqs.(17a) and (19) of
!! 1206.4998.
!!
!! The coefficient functions for the total V and H cross sections have
!! been taken from hep-ph/0207004 and Ellis, Stirling & Webber (with
!! some care about factors of x included in the parton luminosity).
module coefficient_functions_VH
  use types; use consts_dp; use convolution_communicator; use convolution
  use dglap_objects
  use qcd
  implicit none
  private

  public :: InitCoeffMatrix
  
  public                        :: InitHiggsCoeffs
  logical,         save, public :: higgs_init_done = .false.
  type(grid_conv), save, public :: higgs_H1gg, higgs_H1qg, higgs_H1qqbar
  public                        :: InitDYCoeffs
  logical,         save, public :: DY_init_done = .false.
  type(grid_conv), save, public :: DY_D1qqbar, DY_D1qg

contains
   !----------------------------------------------------------------------
   function C1qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = CF*(1-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = -pisq/12._dp*CF
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1qq

   !----------------------------------------------------------------------
   function C1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = CF*x
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1gq


   !----------------------------------------------------------------------
   function C1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 2*TR*x*(1-x)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1qg


  !----------------------------------------------------------------------
  function C1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = 0.0_dp
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = -pisq/12._dp*CA
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function C1gg

   
  !======================================================================
  !! Initialise a coefficient function matrix, repeating exactly
  !! what's done for a LO unpolarised splitting matrix, with the nf
  !! value that is current from the qcd module.
  subroutine InitCoeffMatrix(grid, C)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: C

    !-- info to describe the splitting function
    !C%loops  = 1
    C%nf_int = nf_int

    call cobj_InitSplitLinks(C)

    call InitGridConv(grid, C%gg, C1gg)
    call InitGridConv(grid, C%qq, C1qq)
    call InitGridConv(grid, C%gq, C1gq)
    call InitGridConv(grid, C%qg, C1qg)

    !-- now fix up pieces so that they can be directly used as a matrix
    call Multiply(C%qg, 2*nf)

    !-- PqqV +- PqqbarV
    call InitGridConv(C%NS_plus,  C%qq)
    call InitGridConv(C%NS_minus, C%qq)

    !-- PNSminus + nf * (PqqS - PqqbarS)
    call InitGridConv(C%NS_V, C%NS_minus)
  end subroutine InitCoeffMatrix




  !----------------------------------------------------------------------
  ! form taken from hep-ph/0207004, but dates back much earlier
  ! (e.g. Djouadi et al).
  function H1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      !! expression in Anastasiou & Melnikov hep-ph/0207004
      res = 12.0_dp*log(1-x)*(one/(1-x) - x*(2 - x + x**2))&
           & - 6*(x**2 + 1 - x)**2/(1-x)*log(x)&
           & - 11.0_dp/2.0_dp*(1-x)**3
      !! expression in Dawson, NPB359 (91) 283 -- gives identical results
      ! res = - 11.0_dp/2.0_dp*(1-x)**3 &
      !      & + 6*(1 + x**4 + (1-x)**4)*log(1-x)/(1-x)&
      !      & -6*(x**2/(1-x)+(1-x) + x**2*(1-x))*log(x)
      !res = 12.0_dp*log(1-x)*(one/(1-x) - x*(2 - x*(1-x)))&
      !     & - 11.0_dp/2.0_dp*(1-x)**3
     ! write(0,*) 'REAL', res
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      ! *** NOTE the extra factor of x here and its absence down below.
      ! it's as if the real expression were really "...."/x, to be multiplied
      ! by x below; the 1/x is outside the plus function but the x is not.
      res = res - 12.0_dp*log(1-x)/(1-x)*x
     ! write(0,*) 'REALVIRT', res, 12.0_dp*log(1-x)*(one/(1-x) - x*(2 - x + x**2))&
     !      & - 6*(x**2 + 1 - x)**2/(1-x)*log(x)&
     !      & - 11.0_dp/2.0_dp*(1-x)**3, - 12.0_dp*log(1-x)/(1-x)*x
    case(cc_DELTA)
       res = 11.0_dp/2.0_dp + pisq
    end select

    ! factor of x that usually appears in splitting function is not
    ! present here, cf the note above.
    !
    !!if (cc_piece /= cc_DELTA) res = res * x     

    ! it multiplies alphas/pi, so we need an additional factor of two
    ! for our convention of alphas/twopi
    res = res * two
    !write(0,*) cc_piece, x, res, -12.0_dp*log(1-x)/(1-x)*x + 12.0_dp*log(1-x)*(one/(1-x) - x*(2 - x + x**2))&
!           & - 6*(x**2 + 1 - x)**2/(1-x)*log(x)&
!           & - 11.0_dp/2.0_dp*(1-x)**3
  end function H1gg
  
  !----------------------------------------------------------------------
  ! form taken from hep-ph/0207004, but dates back much earlier
  ! (e.g. Djouadi et al).
  function H1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = -two/three*(1+(1-x)**2) * log(x/(1-x)**2) &
           & - one + two*x - x**2/three
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      ! no virtual piece
    case(cc_DELTA)
      res = zero
    end select
    
    ! it multiplies alphas/pi, so we need an additional factor of two
    ! for our convention of alphas/twopi
    res = res * two
  end function H1qg

  !----------------------------------------------------------------------
  ! form taken from hep-ph/0207004, but dates back much earlier
  ! (e.g. Djouadi et al).
  function H1qqbar(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = 32.0_dp/27.0_dp * (1-x)**3
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      ! no virtual piece
    case(cc_DELTA)
      res = zero
    end select
    
    ! it multiplies alphas/pi, so we need an additional factor of two
    ! for our convention of alphas/twopi
    res = res * two
  end function H1qqbar

  !----------------------------------------------------------------------
  ! This function, and those that follow, are the NLO coefficient
  ! functions for the total Z production cross section, taken from
  ! Ellis, Stirling & Webber.
  function D1qqbar(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = CF*(4*(1+x**2)*log(1-x)/(1-x) - 2*(1+x**2)/(1-x)*log(x))
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      res = res - CF*8*log(1-x)/(1-x)
    case(cc_DELTA)
      res = CF*(2*pisq/three - 8)
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function D1qqbar

  !----------------------------------------------------------------------
  function D1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
       res = TR*((x**2+(1-x)**2)*log((1-x)**2/x) + half + 3*x - 3.5_dp*x**2)
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
       res = zero
    end select

    if (cc_piece /= cc_DELTA) res = res * x     
  end function D1qg

  !--------------------------------------
  subroutine InitHiggsCoeffs(grid)
    type(grid_def), intent(in) :: grid
    if (.not. higgs_init_done) then
      call InitGridConv(grid, higgs_H1gg, H1gg)
      call InitGridConv(grid, higgs_H1qg, H1qg)
      call InitGridConv(grid, higgs_H1qqbar, H1qqbar)
      higgs_init_done = .true.
    end if
  end subroutine InitHiggsCoeffs
  
  !--------------------------------------
  subroutine InitDYCoeffs(grid)
    type(grid_def), intent(in) :: grid
    if (.not. DY_init_done) then
      call InitGridConv(grid, DY_D1qg, D1qg)
      call InitGridConv(grid, DY_D1qqbar, D1qqbar)
      DY_init_done = .true.
    end if
  end subroutine InitDYCoeffs
  
end module coefficient_functions_VH
