!========================================================
!--------------------------------------------------------
! Includes a module for the non-inclusive correction
! -------------------------------------------------------
!========================================================
module emsn_tools
  use types; use consts_dp
  use rad_tools; use hoppet_v1
  implicit none
  private

  public :: non_incl, non_incl_lnR, av_lnz_smallR, av_ln2z_smallR, non_incl_aslnR_sq

  real(dp), parameter :: lntwo =&
       & 0.693147180559945309417232121458176568075_dp

contains

  ! contains the analytic expressions for the non-inclusive (R-dependent) 
  ! corrections
  function non_incl(R,part) result(res)
    real(dp), intent(in) :: R
    real(dp) :: correl, indep, res 
    character(len=*) :: part

    ! Taylor series of the correlated emission contribution neglecting terms of O(R**8)
    ! Good agreement with MC below R=2.6
    !
    ! Taken from 1203.5773 A.16 (and independently derived by PFM)
    correl = -(-131._dp+12._dp*pisq+132._dp*lntwo)/72._dp*ca_def*log(R/1.74_dp) & 
         & -(23.0_dp-24.0_dp*lntwo)/72._dp*nf_def*log(R/0.84_dp) &
         & +((1429._dp+3600._dp*pisq+12480._dp*lntwo)*ca_def+&
         & (3071._dp-1680._dp*lntwo)*nf_def)/172800._dp*R**2           &
         & +((-9383279._dp-117600._dp*pisq+1972320._dp*lntwo)*ca_def+&
         & two*(178080._dp*lntwo-168401._dp)*nf_def)/406425600._dp*R**4 &
         &  +((74801417._dp-33384960._dp*lntwo)*ca_def+&
         & (7001023._dp-5322240._dp*lntwo)*nf_def)/97542144000._dp*R**6

    if (R>pi) then
       ! Taken from 1206.4998 Eq.(45)
       indep = -half*A(1)*((pi/6._dp*R**2-1/pi/8._dp*R**4)*&
            & atan(pi/sqrt(R**2-pisq)) &!this avoids Floating point exceptions for R<pi
            & +(R**2/8._dp-pisq/12._dp)*sqrt(R**2-pisq)) 
    else
       ! Taken from 1203.5773 A.18 (and independently derived by PFM)
       indep = -half*A(1)*(pisq*R**2/12._dp-R**4/16._dp)
    end if

    select case(trim(part))
    case('correl')
       res = correl
    case('indep')
       res = indep
    case('all')
       res = correl+indep
    end select
  end function non_incl

  ! - Added by FD -
  ! contains the analytic expressions for the non inclusive log R correction
  function non_incl_lnR(R,R0) result(res)
    real(dp), intent(in) :: R, R0
    real(dp) :: res 
    ! res = (-131._dp+12._dp*pisq+132._dp*lntwo)/144._dp*ca_def*log(1._dp/(R**2)) & 
    !      & +(23.0_dp-24.0_dp*lntwo)/144._dp*nf_def*log(1._dp/(R**2))
    res = -log(R/R0)*( (23._dp - 24._dp*lntwo)*nf_def &
         & + (-131._dp + 12._dp*pisq + 132._dp*lntwo)*ca_def )/72._dp 
  end function non_incl_lnR
    
  ! - Added by FD -
  ! Return the average ln z, from a parametrized fit to the all-order curve
  ! The first three orders in t are exact, higher orders are the result
  ! of the fit, and do not necessarily correspond to the physical expansion.
  function av_lnz_smallR(as2pi,R,R0) result(res)
    real(dp), intent(in) :: as2pi(:), R, R0
    real(dp) :: t(size(as2pi)),b0
    real(dp) :: aslnR2pi(size(as2pi))
    real(dp) :: res(size(as2pi))
    integer :: i
    ! v1: up to alphas^2
    ! aslnR2pi=as2pi*log(1.0_dp/(R*R))
    ! res = -aslnR2pi*((-131._dp+12._dp*pisq+132._dp*lntwo)/72._dp*ca_def & 
    !      & +(23.0_dp-24.0_dp*lntwo)/72._dp*nf_def) &
    !      & ((aslnR2pi**2)/(2._dp))*(0.3698206876035498_dp*cf_def*nf_def*0.5_dp&
    !      & -0.11786051234372803_dp*(nf_def**2)*0.25_dp-0.5892369676620957_dp*ca_def*nf_def*0.5_dp&
    !      & +0.9015682990280285_dp*ca_def**2)
             
    ! v2: up to alphas^5, where alphas^4,5 are fitted to resummed curve
    ! aslnR2pi = as2pi*log(R)
    ! res = (  aslnR2pi*( (23._dp - 24._dp*lntwo)*nf_def &
    !      & + (-131._dp + 12._dp*pisq + 132._dp*lntwo)*ca_def )/36._dp &
    !      & + (aslnR2pi**2)*( -3.60627_dp*(ca_def**2) + 1.17847_dp*ca_def*nf_def &
    !      & - 0.739641_dp*cf_def*nf_def + 0.117861_dp*(nf_def**2) ) &
    !      & + (aslnR2pi**3)*( 8.57984_dp*(ca_def**3) - 4.97749_dp*(ca_def**2)*nf_def &
    !      & + 2.7538_dp*ca_def*cf_def*nf_def + 0.527541_dp*ca_def*(nf_def**2) &
    !      & - 0.361327_dp*nf_def*(cf_def**2) - 0.636779_dp*cf_def*(nf_def**2) &
    !      & + 0.0785739_dp*(nf_def**3) ) &
    !      & + (aslnR2pi**4)*(-2.02532_dp) + (aslnR2pi**5)*(1.92405_dp)  )
    ! Orders 4 and 5 assume nf=5, cf=4/3, ca=3

    ! v3: using a fit of the resummed curve as a function of t
    !     with the complete expression for t^1,2,3
    ! for ln R resummation, calculate b0 (note different norm from beta0) and t
    b0 = (11._dp*ca_def-2._dp*nf_def)/6._dp
    t = log(one/(one + two*as2pi*b0*log(R/R0)))/b0
    res = ( t*( (-23._dp + 24._dp*lntwo)*nf_def &
         & + (131._dp -12._dp*pisq - 132._dp*lntwo)*ca_def)/72._dp &
         & + (t**2)*( 0.103336_dp*(ca_def**2) + 0.192938_dp*ca_def*nf_def &
         & - 0.18491_dp*cf_def*nf_def + 0.0147326_dp*(nf_def**2) ) &
         & + (t**3)*( -0.0337133_dp*(ca_def**3) - 0.0446767_dp*(ca_def**2)*nf_def &
         & - 0.00522325_dp*ca_def*cf_def*nf_def + 0.0451658_dp*(cf_def**2)*nf_def &
         & - 0.0240506_dp*ca_def*(nf_def**2) + 0.0179606*cf_def*(nf_def**2) &
         & - 0.00163696_dp*(nf_def**3) ))

    ! add higher order coefficients fitted to all order curve, for cf=4/3, ca=3
    if (nf_def.eq.5) then
       ! nf = 5 case
       res = res + (t**4)*(5.56441_dp) + (t**5)*(-3.98790_dp) + (t**6)*(-1.70399_dp) &
            &    + (t**7)*(4.47420_dp) + (t**8)*(-1.91022_dp)
       ! <ln z>_g =  -3.73075528158 t        +  5.91934175      t^2/2     +  -24.2005852431  t^3/6
       !          +  133.545867568  t^4/24   +  -478.547944772  t^5/120   +  -1226.86986449  t^6/720
       !          +  22549.9916024  t^7/5040 +  -77020.0849334  t^8/40320


    elseif (nf_def.eq.4) then
       ! nf = 4 case
       res = res + (t**4)*(4.19559_dp) + (t**5)*(-2.93417_dp) + (t**6)*(-1.19228_dp) &
            &    + (t**7)*(3.13888_dp) + (t**8)*(-1.32930_dp)       
       ! <ln z>_g =  -3.64235989732 t        +  4.989622        t^2/2     +  -18.9422908444  t^3/6
       !          +  100.693874578  t^4/24   +  -352.099959729  t^5/120   +  -858.44000417   t^6/720
       !          +  15819.965337   t^7/5040 +  -53597.497853   t^8/40320

    endif
    
    ! MD+GPS alternative parametrization that works a little better
    ! at high t (but still not perfectly; and probably worse at low t).
    !     res = ( t*( (-23._dp + 24._dp*lntwo)*nf_def &
    !          & + (131._dp -12._dp*pisq - 132._dp*lntwo)*ca_def)/72._dp &
    !          & + (t**2)*( 0.103336_dp*(ca_def**2) + 0.192938_dp*ca_def*nf_def &
    !          & - 0.18491_dp*cf_def*nf_def + 0.0147326_dp*(nf_def**2) ) &
    !          & + (t**3)*( -0.0337133_dp*(ca_def**3) - 0.0446767_dp*(ca_def**2)*nf_def &
    !          & - 0.00522325_dp*ca_def*cf_def*nf_def + 0.0451658_dp*(cf_def**2)*nf_def &
    !          & - 0.0240506_dp*ca_def*(nf_def**2) + 0.0179606*cf_def*(nf_def**2) &
    !          & - 0.00163696_dp*(nf_def**3) ) &
    !          & + (t**4)*3.85003_dp)/(one + (t**5)*(-0.548756_dp) )
    
    
  end function av_lnz_smallR


  
  ! - Added by FD -
  ! term in (as lnR**2)**2 for fixed-order expansion
  function non_incl_aslnR_sq(R, R0) result(res)
    real(dp), intent(in) :: R, R0
    real(dp) :: res 
    ! we have
    !    h_31 = 16 * C * res
    res = log(R/R0)**2 * ( 1.803135_dp*(ca_def**2) - 0.5892350_dp*ca_def*nf_def &
         & + 0.3698205_dp*cf_def*nf_def - 0.0589305_dp*(nf_def**2) )
  end function non_incl_aslnR_sq


  ! - Added by MD & GPS -
  ! contains the average ln^2 z (parametrised up to a third power of t; probably
  ! quite a rough parametrization -- at best a few percent)
  function av_ln2z_smallR(as2pi,R,R0) result(res)
    real(dp), intent(in) :: as2pi(:), R, R0
    real(dp) :: t(size(as2pi)),b0
    real(dp) :: aslnR2pi(size(as2pi))
    real(dp) :: res(size(as2pi))
    integer :: i

    b0 = (11._dp*ca_def-2._dp*nf_def)/6._dp
    t = log(one/(one + two*as2pi*b0*log(R/R0)))/b0
    res = 1.9432_dp*t + 6.1125_dp*t**2 - 2.03343_dp*t**3

    !res = 1.9432_dp * two*as2pi * log(R/R0)
    
  end function av_ln2z_smallR

end module emsn_tools
