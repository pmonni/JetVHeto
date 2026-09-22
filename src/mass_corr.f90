! This module contains the ingredients to compute the top/bottom-quark
! mass corrections in the hard virtual coefficient to NLO
module mass_corr
  use types; use consts_dp
  use ew_parameters
  use hoppet_v1
  implicit none
  private

  public :: c1h, norm_mass

contains
  !======================================================================
  ! running quark mass as defined in Spira, Djouadi, Graudenz & Zerwas
  ! hep-ph/9504378
  ! Different from usual MSbar mass
  function mq(mpole,mu) result(res)
    real(dp), intent(in) :: mpole, mu
    real(dp)             :: res
    !--------------------------------------------- 
    res = mpole*(RunningCoupling(mu)/RunningCoupling(mpole))**(one/pi/beta0)
  end function mq

  !======================================================================
  ! correction factor to account for mass effects in the Born 
  ! cross section for Higgs production
  ! expressions from Harlander & Kant hep-ph/0509189
  ! norm_mass is the ratio of the "finite-mass" cross section
  ! to the "heavy-quark" result
  ! Agreement with hnnlo-v2.0 if muq=mpole (q=t,b)
  ! mut is the scale at which the top-quark mass is defined
  ! mub is the scale at which the bottom-quark mass is defined
  function norm_mass(M, mut, mub, case) result(res)
    character(len=*), intent(in) :: case
    real(dp)        , intent(in) :: M, mut, mub
    !---------------------------------------------
    real(dp) :: res, taub, taut

    ! Eq. (2.2) 
    taub = M**2/four/mq(mb,mub)**2
    taut = M**2/four/mq(mt,mut)**2
    select case(trim(case))
    case('t')
       res = abs(f0(taut))**2
    case('b')
       res = abs(f0(taub))**2
    case('t+b')
       res = abs(f0(taut) + f0(taub))**2
    case default
       call wae_error('heavy-quark mass corrections:  unrecognized case: '//case)
    end select
  end function norm_mass

  !======================================================================
  ! two-loop virtual corrections with heavy quarks running in the loop
  ! expressions from Harlander & Kant hep-ph/0509189
  ! Eqs. (3.5,3.6)
  ! mut is the scale at which the top-quark mass is defined
  ! mub is the scale at which the bottom-quark mass is defined
  function c1h(M,mut,mub,case) result(res)
    character(len=*), intent(in) :: case
    real(dp)        , intent(in) :: M, mut, mub
    !---------------------------------------------
    real(dp) :: res, taub, taut

    taub = M**2/four/mq(mb,mub)**2
    taut = M**2/four/mq(mt,mut)**2
    select case(trim(case))
    case('t')
       res = real((f0b1(theta(taut)) + two*f0c2(taut)*log(mut**2/mq(mt,mut)**2))/f0(taut))
    case('b')
       res = real((f0b1(theta(taub)) + two*f0c2(taub)*log(mub**2/mq(mb,mub)**2))/f0(taub))
    case('t+b')
       res = real((f0b1(theta(taut)) + two*f0c2(taut)*log(mut**2/mq(mt,mut)**2) &
            &    + f0b1(theta(taub)) + two*f0c2(taub)*log(mub**2/mq(mb,mub)**2))/(f0(taut)+f0(taub)))
    case default
       call wae_error('heavy-quark mass corrections:  unrecognized case: '//case)
    end select

    res = res + pisq

  end function c1h

  !======================================================================
  ! Eq. (2.4) 
  function f(tau) result(res)
    real(dp), intent(in) :: tau
    complex(dp) :: res
    
    if(tau < one) then
       res = asin(sqrt(tau))**2
    else
       res = -one/four*(log((one+sqrt(one-one/tau))/(one-sqrt(one-one/tau)))-dcmplx(zero,one)*pi)**2
    endif

  end function f
  !======================================================================
  ! Derivative of Eq. (2.4) 
  function f_p(tau) result(res)
    real(dp), intent(in) :: tau
    complex(dp) :: res
    
    if (tau < one) then
       res = asin(sqrt(tau))/sqrt(tau*(one-tau))
    else 
       res = -half*(log((one+sqrt(one-one/tau))/(one-sqrt(one-one/tau)))-dcmplx(zero,one)*pi)&
            & *(one-sqrt(one-one/tau))/(one+sqrt(one-one/tau))&
            & *one/((sqrt(tau-one)-sqrt(tau))**2*sqrt(tau*(tau-one)))
    endif

  end function f_p
  !======================================================================
  ! Eq. (2.10) 
  function theta(tau) result(res) 
    real(dp), intent(in) :: tau
    complex(dp) :: res

!    res = (sqrt(one-one/tau)-one)/(sqrt(one-one/tau)+one)
    if (tau < one) then
       res =dcmplx(one-two*tau,two*sqrt(tau*(one-tau)))
    else
       res = dcmplx((sqrt(tau-one)-sqrt(tau))/(sqrt(tau-one)+sqrt(tau)), zero)
    end if
  end function theta


  !------------------------------
!!$  function f0b1(tau) result(res)
!!$    real(dp), intent(in) :: tau
!!$    complex(dp) :: res
!!$
!!$    if (tau < one) then
!!$       res = f0b1_chaplin(theta(tau))
!!$    else
!!$       res = f0b1_hplog(real(theta(tau)))
!!$    end if
!!$
!!$  end function f0b1

  !======================================================================
  ! Eq. (3.7)
!!$  function f0b1_hplog(t) result(res)
!!$    real(dp), intent(in) :: t
!!$    complex(dp)          :: res
!!$    !------------------------------------------
!!$    integer, parameter :: n1=-1
!!$    integer, parameter :: n2= 1
!!$    integer, parameter :: nw= 4
!!$    real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
!!$         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
!!$    real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
!!$         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
!!$    complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
!!$         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
!!$
!!$    call hplog(t,nw,Hc1,Hc2,Hc3,Hc4, &
!!$         &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
!!$
!!$    res  =  (-94d0*t)/(1 - t)**2 + (2d0*pi**4*t*(1 + t - 3d0*t**2 -&
!!$         &     3d0*t**3))/&
!!$         &     (5d0*(1 - t)**5) - (24d0*t*(1 + t)*Hc1(0))/(1 - t)**3 + &
!!$         &     (3d0*t*(3d0 + 22d0*t + 3d0*t**2)*Hc1(0)**2)/(2d0*(1 - t)**4)+ &
!!$         &     (pi**2*t*(1 + t - 17d0*t**2 - 17d0*t**3)*Hc1(0)**2)/&
!!$         &     (6d0*(1 - t)**5) + &
!!$         &     (t*(6d0 + 59d0*t + 58d0*t**2 + 33d0*t**3)*Hc1(0)**3)/&
!!$         &     (3d0*(1 - t)**5)+ &
!!$         &     (t*(5d0 + 5d0*t - 13d0*t**2 - 13d0*t**3)*Hc1(0)**4)&
!!$         &     /(24d0*(1 - t)**5)- &
!!$         &     (t*(47d0 + 66d0*t + 47d0*t**2)*Hc1(0)**2*Hc1(1))/(1 - t)**4 + &
!!$         &     (16d0*t*(1 + t + t**2 + t**3)*Hc1(0)**2*Hc2(0,-1))/(1 - t)**5&
!!$         &     + (2*t*(51d0 + 74d0*t + 51d0*t**2)*Hc1(0)*Hc2(0,1))/&
!!$         &     (1 - t)**4 - &
!!$         &     (2d0*t*(5d0 + 5d0*t + 23d0*t**2 + 23d0*t**3)*&
!!$         &     Hc1(0)**2*Hc2(0,1))/(1 - t)**5 - &
!!$         &     (2d0*t*(55d0 + 82d0*t + 55d0*t**2)*Hc3(0,0,1))/(1 - t)**4 + &
!!$         &     (4d0*t*(1 + t)*(23d0 + 41d0*t**2)*Hc1(0)*&
!!$         &     (-Hc3(0,0,-1) + Hc3(0,0,1)))/&
!!$         &     (1 - t)**5 + (36d0*t*(5d0 + 5d0*t + 11d0*t**2 + 11d0*t**3)*&
!!$         &     Hc4(0,0,0,-1))/&
!!$         &     (1 - t)**5 - (36d0*t*(5d0 + 5d0*t + 7d0*t**2 + 7d0*t**3)*&
!!$         &     Hc4(0,0,0,1))/&
!!$         &     (1 - t)**5 + (2d0*t*(31d0 + 34d0*t + 31d0*t**2)*zeta3)/&
!!$         &     (1 - t)**4 + &
!!$         &     (2d0*t*(11d0 + 11d0*t - 43d0*t**2 - 43d0*t**3)*Hc1(0)*zeta3)/&
!!$         &     (1 - t)**5 + &
!!$         &     (t*(1 + t)**2*((-4d0*pi**2*Hc1(0))/3d0 + &
!!$         &     6d0*pi**2*Hc1(0)*Hc1(1) - &
!!$         &     6d0*Hc1(0)**3*Hc1(1) - 32d0*Hc1(0)*Hc2(0,-1) -&
!!$         &     6d0*pi**2*Hc2(0,1) + &
!!$         &     64d0*Hc3(0,0,-1) + 72d0*Hc4(1,0,-1,0) + 108d0*Hc1(1)*zeta3))/&
!!$         &     (1 - t)**4
!!$
!!$! GZ check     
!!$!    res2 = (t*(11280*(-1 + t)**3 + 48*Pi**4*(1 + t - 3*t**2 - 3*t**3) -                     &
!!$!     &      240*(-1 + t)*(31 + 34*t + 31*t**2)*zeta3 +                                      &
!!$!     &      240*(-1 + t)*(55 + 82*t + 55*t**2)*Hc3(0,0,1) +                                 &
!!$!     &      4320*(5 + 5*t + 11*t**2 + 11*t**3)*Hc4(0,0,0,-1) -                              &
!!$!     &      4320*(5 + 5*t + 7*t**2 + 7*t**3)*Hc4(0,0,0,1) -                                 &
!!$!     &      2880*(-1 + t)**2*(1 + t)*Log(t) +                                               &
!!$!     &      240*(11 + 11*t - 43*t**2 - 43*t**3)*zeta3*Log(t) -                              &
!!$!     &      240*(-1 + t)*(51 + 74*t + 51*t**2)*Hc2(0,1)*Log(t) +                            &
!!$!     &      480*(1 + t)*(23 + 41*t**2)*(-Hc3(0,0,-1) + Hc3(0,0,1))*Log(t) -                 &
!!$!     &      180*(-1 + t)*(3 + 22*t + 3*t**2)*Log(t)**2 +                                    &
!!$!     &      20*Pi**2*(1 + t - 17*t**2 - 17*t**3)*Log(t)**2 +                                &
!!$!     &      120*(-1 + t)*(47 + 66*t + 47*t**2)*Hc1(1)*Log(t)**2 +                           &
!!$!     &      1920*(1 + t + t**2 + t**3)*Hc2(0,-1)*Log(t)**2 -                                &
!!$!     &      240*(5 + 5*t + 23*t**2 + 23*t**3)*Hc2(0,1)*Log(t)**2 +                          &
!!$!     &      40*(6 + 59*t + 58*t**2 + 33*t**3)*Log(t)**3 +                                   &
!!$!     &      5*(5 + 5*t - 13*t**2 - 13*t**3)*Log(t)**4 -                                     &
!!$!     &      80*(-1 + t)*(1 + t)**2*(3*                                                      &
!!$!     &          (54*zeta3*Hc1(1) - 3*Pi**2*Hc2(0,1) + 32*Hc3(0,0,-1) + 36*Hc4(1,0,-1,0)) +  &
!!$!     &         (Pi**2*(-2 + 9*Hc1(1)) - 48*Hc2(0,-1))*Log(t) - 9*Hc1(1)*Log(t)**3)))/       &
!!$!     &  (120.*(1 - t)**5)
!!$!   write(*,*) 'res',res-res2,res,res2
!!$
!!$  end function f0b1_hplog

  function f0b1(t) result(res)
    complex(dp), intent(in) :: t
    complex(dp)          :: res
    complex(dp) :: HPL1, HPL2, HPL3, HPL4


    res  =  (-94d0*t)/(1 - t)**2 + (2d0*pi**4*t*(1 + t - 3d0*t**2 -&
         &     3d0*t**3))/&
         &     (5d0*(1 - t)**5) - (24d0*t*(1 + t)*HPL1(0,t))/(1 - t)**3 + &
         &     (3d0*t*(3d0 + 22d0*t + 3d0*t**2)*HPL1(0,t)**2)/(2d0*(1 - t)**4)+ &
         &     (pi**2*t*(1 + t - 17d0*t**2 - 17d0*t**3)*HPL1(0,t)**2)/&
         &     (6d0*(1 - t)**5) + &
         &     (t*(6d0 + 59d0*t + 58d0*t**2 + 33d0*t**3)*HPL1(0,t)**3)/&
         &     (3d0*(1 - t)**5)+ &
         &     (t*(5d0 + 5d0*t - 13d0*t**2 - 13d0*t**3)*HPL1(0,t)**4)&
         &     /(24d0*(1 - t)**5) - &
         &     (t*(47d0 + 66d0*t + 47d0*t**2)*HPL1(0,t)**2*HPL1(1,t))/(1 - t)**4 + &
         &     (16d0*t*(1 + t + t**2 + t**3)*HPL1(0,t)**2*HPL2(0,-1,t))/(1 - t)**5 &
         &     + (2*t*(51d0 + 74d0*t + 51d0*t**2)*HPL1(0,t)*HPL2(0,1,t))/&
         &     (1 - t)**4 - &
         &     (2d0*t*(5d0 + 5d0*t + 23d0*t**2 + 23d0*t**3)*&
         &     HPL1(0,t)**2*HPL2(0,1,t))/(1 - t)**5 - &
         &     (2d0*t*(55d0 + 82d0*t + 55d0*t**2)*HPL3(0,0,1,t))/(1 - t)**4 + &
         &     (4d0*t*(1 + t)*(23d0 + 41d0*t**2)*HPL1(0,t)*&
         &     (-HPL3(0,0,-1,t) + HPL3(0,0,1,t)))/&
         &     (1 - t)**5 + (36d0*t*(5d0 + 5d0*t + 11d0*t**2 + 11d0*t**3)*&
         &     HPL4(0,0,0,-1,t))/&
         &     (1 - t)**5 - (36d0*t*(5d0 + 5d0*t + 7d0*t**2 + 7d0*t**3)*&
         &     HPL4(0,0,0,1,t))/&
         &     (1 - t)**5 + (2d0*t*(31d0 + 34d0*t + 31d0*t**2)*zeta3)/&
         &     (1 - t)**4 + &
         &     (2d0*t*(11d0 + 11d0*t - 43d0*t**2 - 43d0*t**3)*HPL1(0,t)*zeta3)/&
         &     (1 - t)**5 + &
         &     (t*(1 + t)**2*((-4d0*pi**2*HPL1(0,t))/3d0 + &
         &     6d0*pi**2*HPL1(0,t)*HPL1(1,t) - &
         &     6d0*HPL1(0,t)**3*HPL1(1,t) - 32d0*HPL1(0,t)*HPL2(0,-1,t) -&
         &     6d0*pi**2*HPL2(0,1,t) + &
         &     64d0*HPL3(0,0,-1,t) + 72d0*HPL4(1,0,-1,0,t) + 108d0*HPL1(1,t)*zeta3))/&
         &     (1 - t)**4

  end function f0b1

  !======================================================================
  ! Eq. (2.7) 
  function f0c2(t) result(res)
    real(dp), intent(in) :: t
    complex(dp)          :: res
    !------------------------------------------
       
    res = three/t**2*(t+(t-two)*f(t)&
         & - (t-one)*t*f_p(t))
  end function f0c2
  !======================================================================
  ! Eq. (2.3) 
  function f0(t) result(res)
    real(dp), intent(in) :: t
    complex(dp)          :: res
    !------------------------------------------
    
    res = half*three/t**2*(t+(t-one)*f(t))
   end function f0

end module mass_corr
