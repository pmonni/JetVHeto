module resummation
  use types; use consts_dp
  use rad_tools
  use pdfs_tools
  use warnings_and_errors
  use emsn_tools
  use ew_parameters; use mass_corr
  use special_functions 
  implicit none

  private
  public :: resummed_sigma

contains
  !======================================================================
  function resummed_sigma(pt, cs, order,dlumi_lumi) result(sigma)
    real(dp),                  intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order
    real(dp), intent(out), optional       :: dlumi_lumi(size(pt))
    real(dp)  :: sigma(size(pt))
    !------------------------------
    real(dp) :: L_tilde(size(pt)), lambda(size(pt)), lNLL(size(pt)), &
         &lNNLL(size(pt)), dlNLL(size(pt))
    real(dp) :: normalisation
    real(dp) :: rp(size(pt)),rs(size(pt)),resum_fact(size(pt)),drp(size(pt))
    real(dp) :: av_lnz(size(pt)), av_ln2z(size(pt)), as2pi_pt(size(pt)), non_incl_largeR(size(pt))
    integer :: i 
    real(dp) :: tmp(size(pt))
    

    L_tilde = Ltilde(pt/cs%Q, cs%p)
    lambda  = get_lambda(L_tilde, cs)

    sigma = zero
    where (lambda < half)  sigma = exp(Rad(L_tilde, cs, order))

    if (present(dlumi_lumi))  dlumi_lumi = zero    

    if (order == order_LL) then
       sigma = sigma * lumi_LL(cs)
    else if (order == order_NLL) then
       sigma = sigma * lumi_NLL(exp(-L_tilde)*cs%muF, cs)

       if (cs%observable == 'ptB') then 
          rp=Rad_p(lambda)
          do i=1,size(lambda)
             sigma(i) = sigma(i) * &
                  & exp(-eulergamma*rp(i))*gammafull(one-rp(i)/two)&
                  & /gammafull(one+rp(i)/two) 
          enddo
       endif

    else

       if (order /= order_NNLL) call wae_error("expected order_NNLL, found", intval=order)

       ! include the hard correction and the NNLL luminosity
       lNLL  = lumi_NLL(exp(-L_tilde)*cs%muF, cs)
       lNNLL = lumi_NNLL(exp(-L_tilde)*cs%muF, lambda, cs)
       
       if (cs%observable == 'ptj') then 

          if (cs%small_r) then
             ! - Added by FD -
             ! get the all-order result
             as2pi_pt = cs%as2pi/(one-two*lambda)
             ! for ln R resummation, calculate b0 and t
             ! b0 = (11.*ca_def-2.*nf_def)/6.
             ! t = log(1/(1-as2pi_pt*b0*log(1/(cs%jet_radius**2))/(2.*pi)))/b0
             av_lnz = av_lnz_smallR(as2pi_pt,cs%jet_radius,cs%small_r_R0)
             ! subtract lnR term at O(as) from non-inclusive correction
             non_incl_largeR = as2pi_pt*(non_incl(cs%jet_radius,'all') - &
                  & non_incl_lnR(cs%jet_radius,cs%small_r_R0))
             ! include the non-inclusive correction
             ! this should be checked carefully
             sigma = sigma * (exp(-Rad_p(lambda)*av_lnz) + &
                  & Rad_p(lambda)*two*non_incl_largeR)
             ! sigma = sigma * (1 - Rad_p(lambda)*av_lnz + &
             !      & Rad_p(lambda)*two*non_incl_largeR)

             if (cs%small_r_ln2z) then
                ! extra subleading terms: MD & GPS temporary investigations (2015-02-16)
                av_ln2z = av_ln2z_smallR(as2pi_pt,cs%jet_radius,cs%small_r_R0)
                !
                ! fix a normalisation so as not to change results at large pt
                tmp(1:1) = exp(-A(1)*cs%as2pi*two &
                     &          * av_ln2z_smallR((/cs%as2pi/),cs%jet_radius,cs%small_r_R0))
                ! put it together
                sigma = sigma * exp(-A(1)*as2pi_pt*two * av_ln2z) / tmp(1)
                
             end if
          else
             ! - Original code -
             sigma = sigma * (1 + &
                  & non_incl(cs%jet_radius,'all') * Rad_p(lambda)*two*cs%as2pi/(1-two*lambda))
          end if

          if (cs%include_c1_squared) then
             ! Option 1 for luminosity: everything factorised
             if (cs%loop_mass .ne. 'none') then
                sigma = sigma * (1+cs%as2pi*Hm(1))* lNNLL
             else
                sigma = sigma * (1+cs%as2pi*H(1))* lNNLL
             endif
             if (present(dlumi_lumi)) then
                if (cs%loop_mass .ne. 'none') then
                   dlumi_lumi = (1+cs%as2pi*Hm(1))* lNNLL/lNLL-one
                else
                   dlumi_lumi = (1+cs%as2pi*H(1))* lNNLL/lNLL-one
                endif
             end if
          else
             ! Option 2 for luminosity: no spurious O(as^2) terms.
             ! This is the option that we take as our default.
             if (cs%loop_mass .ne. 'none') then
                sigma = sigma * ( lNNLL+ lNLL*cs%as2pi*Hm(1) )
             else
                sigma = sigma * ( lNNLL+ lNLL*cs%as2pi*H(1) )
             endif
             if (present(dlumi_lumi)) then
                if (cs%loop_mass .ne. 'none') then
                   dlumi_lumi = lNNLL/lNLL+ cs%as2pi*Hm(1)-one
                else 
                   dlumi_lumi = lNNLL/lNLL+ cs%as2pi*H(1)-one
                endif
             end if
          end if

       else ! doing ptB 
          dlNLL  = dlumi_NLL(exp(-L_tilde)*cs%muF, cs)
          resum_fact = 1d0 

          ! drp is the NNLL contribution to rp
          drp=Rad_pNNLL(lambda,cs)
          rp=Rad_p(lambda)
          rs=Rad_s(lambda,cs) 
          if (cs%include_c1_squared) then
             ! Option 1 for luminosity: everything factorised
             if (cs%loop_mass .ne. 'none') then
                resum_fact = resum_fact * (1+cs%as2pi*Hm(1))* lNNLL
             else
                resum_fact = resum_fact * (1+cs%as2pi*H(1))* lNNLL
             endif
             if (present(dlumi_lumi)) then
                if (cs%loop_mass .ne. 'none') then
                   dlumi_lumi = (1+cs%as2pi*Hm(1))* lNNLL/lNLL-one
                else
                   dlumi_lumi = (1+cs%as2pi*H(1))* lNNLL/lNLL-one
                endif
             end if
          else
             ! Option 2 for luminosity: no spurious O(as^2) terms.
             ! This is the option that we take as our default.
             if (cs%loop_mass .ne. 'none') then
                resum_fact = resum_fact * ( lNNLL+ lNLL*cs%as2pi*Hm(1) )
             else
                resum_fact = resum_fact * ( lNNLL+ lNLL*cs%as2pi*H(1) )
             endif
             if (present(dlumi_lumi)) then
                if (cs%loop_mass .ne. 'none') then
                   dlumi_lumi = lNNLL/lNLL+ cs%as2pi*Hm(1)-one
                else 
                   dlumi_lumi = lNNLL/lNLL+ cs%as2pi*H(1)-one
                endif
             end if
          end if
    
          do i=1,size(lambda)
             resum_fact(i) = resum_fact(i) * &
                  & exp(-eulergamma*rp(i))*gammafull(one-rp(i)/two)&
                  & /gammafull(one+rp(i)/two)

             resum_fact(i) = resum_fact(i) - rs(i)/two* & 
                  &        ((gammafull(one - rp(i)/two)*                     &
                  &        (-EulerGamma**2 + dpsipg(one - rp(i)/two,0)**2 + &
                  &        dpsipg(one - rp(i)/two,1) +                      &
                  &        (EulerGamma + dpsipg(one + rp(i)/two,0))**2 +    &
                  &      2*(EulerGamma + dpsipg(one - rp(i)/two,0))*        &
                  &       (2*EulerGamma + dpsipg(one + rp(i)/two,0)) -      &
                  &      dpsipg(one + rp(i)/two,1)))/                       &
                  &  (4.*exp(EulerGamma*rp(i))*gammafull(one + rp(i)/two)))*lNLL(i)

             resum_fact(i) = resum_fact(i) -&
                  &  (-(gammafull(1 - rp(i)/two)* &
                  &  (two*EulerGamma + dpsipg(1 - rp(i)/two,0) +     &
                  &   dpsipg(1 + rp(i)/two,0)))/                     & 
                  &  (two*exp(EulerGamma*rp(i))*gammafull(1 + rp(i)/two))) *&
                  & (dlNLL(i) - drp(i)*lNLL(i))
          enddo
          sigma = sigma * resum_fact
       endif
    end if

  end function resummed_sigma


  !! extend dgamma to negative arguments
  function gammafull(x) result (res)
    use types; use consts_dp
    use special_functions
    real(dp)             :: res
    real(dp), intent(in) :: x
    if (x < zero) then
       res = pi/dgamma(one-x)/sin(pi*x)
    else
       res = dgamma(x)
    end if
  end function gammafull


end module resummation
