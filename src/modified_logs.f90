! Shared scale profile and logarithm cutoff for every resummation order.
module modified_logs
  use types
  use rad_tools, only: process_and_parameters, init_proc
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: profile_scale, profile_point
contains
  function profile_scale(pt,M,Q,xM) result(muL)
    real(dp), intent(in) :: pt,M,Q,xM
    real(dp) :: muL,matchscale,t
    if(.not.all(ieee_is_finite([pt,M,Q,xM]))) error stop 'Profile: nonfinite input'
    if(pt<0.or.min(M,Q)<=0.or.xM<0.5_dp.or.xM>1) &
      error stop 'Profile: positive scales and 0.5 <= xM <= 1 required'
    matchscale=xM*M
    if(Q<matchscale/4) error stop 'Profile requires Q >= xM*M/4 to keep the logarithm nonnegative'
    if(pt>=M) then
      muL=M
    else if(pt<matchscale) then
      t=pt/matchscale
      ! RadISH p=2 low-scale polynomial, in dimensionless coordinates.
      muL=pt+(1-t)**3*(Q+(3*Q-matchscale)*t)
    else
      ! Algebraic factorization of RadISH's p=2 quartic bridge.
      ! muL=pt and d(muL)/dpt=1 at both joins, so L and L' vanish there.
      ! This avoids the apparent (M-muM)^4 denominator of the original.
      muL=pt+(pt-M)**2*(pt-matchscale)**2/(2*M*matchscale*(M+matchscale))
    endif
    if(muL<=0.or..not.ieee_is_finite(muL)) error stop 'Profile: invalid profiled scale'
  end function

  subroutine profile_point(cs,pt,local,evalpt)
    type(process_and_parameters), intent(in) :: cs
    real(dp), intent(in) :: pt
    type(process_and_parameters), intent(out) :: local
    real(dp), intent(out) :: evalpt(1)
    local=cs
    local%Q=profile_scale(pt,cs%M,cs%Q,cs%xM)
    local%ln_Q2_M2=2*log(local%Q/local%M)
    local%ln_Q2_muR2=2*log(local%Q/local%muR)
    local%use_new_modlog=.false.
    local%p=-1
    evalpt=min(pt,cs%M)
    ! RadISH new-modlog: L=log(Q_profile/pt)*Theta(xM*M-pt).
    ! Keep Q_profile in the hard/beam coefficients; only switch off L.
    if(pt>cs%xM*cs%M) evalpt=local%Q
    call init_proc(local)
  end subroutine
end module modified_logs
