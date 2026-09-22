module jetveto_n3ll
  use modified_logs, only: profile_scale, profile_point
  use rapidity_n3ll, only: jet_radius_exponent, jet_radius_series
  use types
  use consts_dp
  use qcd
  use rad_tools
  use pdfs_tools, only: lumi_LL
  use luminosity_n3ll
  use process_n3ll
  use radiator_n3ll
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: initialize_jetveto_n3ll, resum_jetveto_n3ll, expand_jetveto_n3ll
  public :: jet_radius_exponent, jet_radius_series
  public :: profile_scale, profile_point
contains
  subroutine initialize_jetveto_n3ll(cs,directory)
    type(process_and_parameters), intent(in) :: cs
    character(len=*), intent(in) :: directory
    real(dp) :: an(4),bn(3),betas(0:3),hard(2),d2,d3
    if(cs%small_r .or. cs%small_r_ln2z .or. cs%include_c1_squared) &
      error stop 'N3LL: small-R resummation and c1-squared options are not supported'
    call coefficients_n3ll(cs,an,bn,betas,hard,d2,d3,cs%jet_algorithm)
    call init_luminosity_n3ll(cs%jet_radius,directory)
  end subroutine

  recursive function resum_jetveto_n3ll(pt,cs,dlumi) result(sigma)
    real(dp), intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    real(dp), optional, intent(out) :: dlumi(size(pt))
    real(dp) :: sigma(size(pt)),L(size(pt)),lam,rad0(1),exponent,pref,born
    real(dp) :: an(4),bn(3),betas(0:3),hard(2),d2,d3
    integer :: i
    type(process_and_parameters) :: local
    real(dp) :: evalpt(1),localdl(1),rfactor
    if(any(pt<=0) .or. .not.all(ieee_is_finite(pt))) error stop 'N3LL requires finite positive veto scales'
    if(cs%use_new_modlog) then
      do i=1,size(pt)
        call profile_point(cs,pt(i),local,evalpt)
        sigma(i:i)=resum_jetveto_n3ll(evalpt,local,localdl)
        if(present(dlumi)) dlumi(i)=localdl(1)
      enddo
      call init_proc(cs)
      return
    endif
    call coefficients_n3ll(cs,an,bn,betas,hard,d2,d3,cs%jet_algorithm)
    born=lumi_LL(cs)
    if(born<=0) error stop 'N3LL: nonpositive Born normalization'
    if(cs%p<0) then
      L=log(cs%Q/pt)
    else
      L=Ltilde(pt/cs%Q,cs%p)
    endif
    sigma=0
    if(present(dlumi)) dlumi=0
    do i=1,size(pt)
      lam=cs%alphas_muR*beta0*L(i)
      if(lam>=half) cycle
      rad0=Rad([L(i)],cs,order_NNLL)
      exponent=rad0(1)+4*cs%as2pi**2*g4_n3ll(lam,an,bn,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2)
      pref=n3ll_lumi(cs%muF*exp(-L(i)),lam,cs,hard(1),hard(2))
      rfactor=exp(jet_radius_exponent(L(i),cs%as2pi,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3, &
        cs%truncate_rapidity_rge_at_n3ll))
      sigma(i)=exp(exponent)*pref*rfactor
      ! anew: keep the radius factor outside the matching prefactor.
      if(present(dlumi)) dlumi(i)=pref/born-1
    end do
    if(.not.all(ieee_is_finite(sigma))) error stop 'N3LL: nonfinite resummed result'
  end function

  recursive function expand_jetveto_n3ll(pt,cs,dlumi) result(v)
    real(dp), intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    real(dp), optional, intent(out) :: dlumi(2,size(pt))
    real(dp) :: v(3,size(pt)),L(size(pt)),an(4),bn(3),betas(0:3),hard(2),d2,d3
    real(dp) :: p(0:3),s(3),f(3),a0,born
    integer :: i
    type(process_and_parameters) :: local
    real(dp) :: evalpt(1),localdl(2,1),prefseries(2)
    if(any(pt<=0).or..not.all(ieee_is_finite(pt))) error stop 'N3LL requires finite positive veto scales'
    if(cs%use_new_modlog) then
      do i=1,size(pt)
        call profile_point(cs,pt(i),local,evalpt)
        v(:,i:i)=expand_jetveto_n3ll(evalpt,local,localdl)
        if(present(dlumi)) dlumi(:,i)=localdl(:,1)
      enddo
      call init_proc(cs)
      return
    endif
    call coefficients_n3ll(cs,an,bn,betas,hard,d2,d3,cs%jet_algorithm)
    if(cs%p<0) then
      L=log(cs%Q/pt)
    else
      L=Ltilde(pt/cs%Q,cs%p)
    endif
    a0=cs%as2pi
    born=lumi_LL(cs)
    do i=1,size(pt)
      p=n3ll_lumi_series(L(i),cs,hard(1),hard(2))
      prefseries=[a0*p(1),a0**2*p(2)]/born
      f=jet_radius_series(L(i),betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3)
      ! F starts at a^2, so its square first enters at a^4.
      p(3)=p(3)+f(2)*p(1)+f(3)*p(0)
      p(2)=p(2)+f(2)*p(0)
      s=radiator_expansion(L(i),an,bn,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2)
      v(1,i)=a0*(p(1)+s(1)*p(0))
      v(2,i)=a0**2*(p(2)+s(1)*p(1)+(s(2)+s(1)**2/2)*p(0))
      v(3,i)=a0**3*(p(3)+s(1)*p(2)+(s(2)+s(1)**2/2)*p(1) &
        +(s(3)+s(1)*s(2)+s(1)**3/6)*p(0))
      if(present(dlumi)) dlumi(:,i)=prefseries
    end do
  end function
end module
