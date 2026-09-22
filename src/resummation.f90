module resummation
  use types; use consts_dp
  use rad_tools
  use pdfs_tools
  use warnings_and_errors
  use emsn_tools
  use ew_parameters; use mass_corr
  use special_functions 
  use modified_logs, only: profile_point
  use jetveto_n3ll, only: resum_jetveto_n3ll
  implicit none

  private
  public :: resummed_sigma

contains
  !======================================================================
  recursive function resummed_sigma(pt, cs, order,dlumi_lumi) result(sigma)
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
    integer :: i 
    real(dp) :: tmp(size(pt))
    type(process_and_parameters) :: local
    real(dp) :: evalpt(1),localdl(1),localL(1),fjet,radlocal(1)
    

    call validate_small_r(cs,order)
    if(order==order_N3LL) then
       sigma=resum_jetveto_n3ll(pt,cs,dlumi_lumi)
       return
    end if
    if(cs%use_new_modlog.or.cs%matching_anew) then
       ! anew: complete evolved luminosity P, with F outside matching.
       do i=1,size(pt)
          if(cs%use_new_modlog) then
             call profile_point(cs,pt(i),local,evalpt)
          else
             local=cs;evalpt=pt(i)
          endif
          local%matching_anew=.false.
          sigma(i:i)=resummed_sigma(evalpt,local,order,localdl)
          localL=Ltilde(evalpt/local%Q,local%p)
          fjet=0
          if(order==order_NNLL) then
             tmp(1:1)=get_lambda(localL,local)
             radlocal=Rad_p(tmp(1:1))
             fjet=non_incl(local%jet_radius,'all')*radlocal(1)*two*local%as2pi/(1-two*tmp(1))
             if(local%small_r) then
                ! The recursive result already contains the full small-R
                ! factor. Remove it only when extracting P for anew.
                radlocal=small_r_factor(tmp(1:1),local)
                fjet=log(radlocal(1))
             else
                ! exp(F2) and 1+F2 have identical expansions through a^3.
                sigma(i)=sigma(i)*exp(fjet)/(1+fjet)
             endif
          endif
          radlocal=Rad(localL,local,order)
          if(present(dlumi_lumi)) dlumi_lumi(i)=sigma(i)/exp(radlocal(1)+fjet)/lumi_LL(local)-1
       enddo
       call init_proc(cs)
       return
    endif
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
             sigma = sigma * small_r_factor(lambda,cs)
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
