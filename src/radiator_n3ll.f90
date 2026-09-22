! N3LL Sudakov g4, imported from RadISH 713d2c015, part0 only.
! Caller supplies A,B in alpha_s/(2*pi) and beta in d alpha_s/d ln(mu^2).
! This includes the RSV zeta3 conversion present in the reference radiator.
module radiator_n3ll
  use types
  use consts_dp
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: g4_n3ll, radiator_expansion
contains
  function g4_n3ll(lambda,A,B,beta,lqm,lqr) result(res)
    real(dp), intent(in) :: lambda,A(4),B(3),beta(0:3),lqm,lqr
    real(dp) :: res,beta0,beta1,beta2,beta3
    if (.not.ieee_is_finite(lambda) .or. lambda<0 .or. lambda>=half) &
      error stop 'N3LL radiator: lambda outside [0,1/2)'
    if (.not.all(ieee_is_finite(A)) .or. .not.all(ieee_is_finite(B)) .or. &
        .not.all(ieee_is_finite(beta)) .or. beta(0)<=0 .or. &
        .not.ieee_is_finite(lqm) .or. .not.ieee_is_finite(lqr)) &
      error stop 'N3LL radiator: invalid coefficient or scale logarithm'
    res=zero
    if(lambda==zero) return
    beta0=beta(0); beta1=beta(1); beta2=beta(2); beta3=beta(3)
    res=(&
               &   2*lambda*(4*beta0**2*pi*(&
               &      2*beta0**4*lqr**3*lambda*(-3 + 2*lambda)*pi**2*A(1)&
               &    - 8*beta0**4*(3 - 6*lambda + 4*lambda**2)*pi**2*zeta3*A(1)&
               &    + 3*beta0**2*lqr**2*pi*(&
               &         beta1*pi*A(1) + beta0*(3 - 2*lambda)*lambda*A(2)&
               &       + 2*beta0**2*(-one+lambda)*(-one+two*lambda)*pi*B(1)   )&
               &    - (3*lqr*(&
               &         2*beta0*beta1*(-one+lambda)*(-one+two*lambda)*pi*A(2)&
               &       + lambda*(-4*(beta1**2-beta0*beta2)*(1 + 2*lambda)*pi**2*A(1) + beta0**2*(3 - 2*lambda)*A(3))&
               &       + 4*beta0**3*(-one+lambda)*(-one+two*lambda)*pi*B(2)   )&
               &      )/2.    )&
               & + 3*beta0*pi*(&
               &      beta0*beta1*A(3) - 2*pi*(2*beta0*beta3*pi*A(1)&
               &    + beta1*(-2*beta2*pi*A(1) + beta1*A(2) + 2*beta0**2*B(2))  )&
               & + 2*beta0**3*B(3)) + 3*beta0*lambda*(&
               &      beta0**2*A(4) + pi*(&
               &        -5*beta0*beta1*A(3) + 2*pi*(3*beta1**2*A(2) + 2*pi*(&
               &            -5*beta1*beta2*A(1) + 5*beta0*beta3*A(1) + 2*beta0*(-beta1**2 + beta0*beta2)*B(1)  )&
               &        + 6*beta0**2*beta1*B(2)  ) - 6*beta0**3*B(3)  )  )&
               & + 2*lambda**2*( -(beta0**3*A(4))&
               &    + pi*(5*beta0**2*beta1*A(3) + 2*pi*(beta0*(-11*beta1**2 + 8*beta0*beta2)*A(2)&
               & + 2*pi*((2*beta1**3 + 5*beta0*beta1*beta2 - 7*beta0**2*beta3)*A(1)&
               & + 6*beta0**2*(beta1**2 - beta0*beta2)*B(1)) - 6*beta0**3*beta1*B(2)) + 6*beta0**4*B(3))))&
               & - 3*pi*(&
               &   beta0**2*beta1*(-1 + 6*lambda)*A(3)&
               & + 4*beta0**2*beta1*lqr*pi*(&
               &   beta0**2*lqr*(-1 + 6*lambda)*pi*A(1) + (beta0 - 6*beta0*lambda)*A(2)&
               & + 2*pi*(-2*beta1*lambda*A(1) + beta0**2*(-one+two*lambda)*B(1)))&
               & + 2*pi*(2*(-(beta0**2*beta3*(-one+two*lambda)**3) - 4*beta1**3*lambda**2*(1 + 2*lambda)&
               & + beta0*beta1*beta2*(-one+two*lambda*(3 - 4*lambda + 8*lambda**2)))*pi*A(1) &
               & + beta0*beta1**2*(one-two*lambda)*A(2) + 2*beta0**3*beta1*(one-two*lambda)*B(2)))*Log(one-two*lambda)&
               & + 6*beta1**2*pi**2*(&
               &   2*beta0**2*lqr*(1 - 6*lambda)*pi*A(1) + beta0*(-1 + 6*lambda)*A(2)&
               & + 2*pi*(2*beta1*lambda*A(1) + beta0**2*(one-two*lambda)*B(1)))*Log(one-two*lambda)**2&
               & + 4*beta1**3*(1 - 6*lambda)*pi**3*A(1)*Log(one-two*lambda)**3&
               & - 12*beta0**2*lqm*(-one+two*lambda)*pi*(lambda*(&
               &   4*beta0**3*lqr*(-one+lambda)*pi*(beta0*lqr*pi*A(1) - A(2))&
               & + 2*pi*(2*(beta1**2 - beta0*beta2)*lambda*pi*A(1) - beta0*beta1*(-one+lambda)*A(2))&
               & + beta0**2*(-one+lambda)*A(3))&
               & + beta0*beta1*pi*(-2*beta0*lqr*pi*A(1) + A(2))*Log(one-two*lambda)&
               & - beta1**2*pi**2*A(1)*Log(one-two*lambda)**2))/(48._dp*beta0**5*(-one+two*lambda)**3*pi**2)
  end function

  ! Independent fixed-order integration of minus the radiator:
  ! -int_0^(2L) dt [ A(a(k))*(t-lqm) + B(a(k)) ],
  ! with k^2=Q^2 exp(-t), plus the order-a^3 RSV conversion.
  ! Returns coefficients of a,a^2,a^3, where a=alpha_s(muR)/(2*pi).
  function radiator_expansion(L,A,B,beta,lqm,lqr) result(c)
    real(dp), intent(in) :: L,A(4),B(3),beta(0:3),lqm,lqr
    real(dp) :: c(3),t,u,b0,b1,running2,running3,f(3),weight
    integer :: i
    c=zero
    b0=twopi*beta(0); b1=twopi**2*beta(1)
    ! Simpson quadrature is exact for the cubic polynomial integrands.
    do i=0,2
      t=L*i; u=lqr-t
      weight=one; if(i==1) weight=four
      running2=-b0*u
      running3=b0**2*u**2-b1*u
      f(1)=A(1)*(t-lqm)+B(1)
      f(2)=(A(2)+running2*A(1))*(t-lqm)+B(2)+running2*B(1)
      f(3)=(A(3)+2*running2*A(2)+running3*A(1))*(t-lqm) &
        +B(3)+2*running2*B(2)+running3*B(1)
      c=c-L/three*weight*f
    end do
    c(3)=c(3)+32*pi**2*beta(0)**2*A(1)*zeta3*L
  end function
end module
