module matching
  use types; use consts_dp
  use warnings_and_errors
  use pdfs_tools; use rad_tools
  implicit none
  
  private
  public :: match, fixed_order_schemes
  integer, public, parameter :: fixed_order_NLO = 1, fixed_order_NNLO = 2, fixed_order_N3LO = 3

contains

  !----------------------------------------------------------------------
  function match(matching_scheme, &
       &            sigma, pt, cs, resummation, resum_expanded, &
       &            sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    character(len=*), intent(in) :: matching_scheme
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    type(process_and_parameters), intent(in) :: cs
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))

    select case(trim(matching_scheme))
    case("a")
       res = matching_a(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("b")
       res = matching_b(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("c")
       res = matching_c(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("cunx")
       res = matching_c_unexp(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("d")
       res = matching_d(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("moda")
       res = matching_moda(sigma, pt(:), cs, resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("modRa")
       res = matching_modRa(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case("R")
       res = matching_R(sigma, pt(:), resummation, resum_expanded, sigmabar) 
    case("loga")
       res = matching_loga(sigma, pt(:), resummation, resum_expanded, &
            &sigmabar(:,:),dlumi_lumi,dlumi_lumi_expanded)
    case default
       call wae_error("matching: unrecognised matching scheme "//matching_scheme)
    end select
  end function match
  
! --- new version
  function matching_a(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    matching_factor = sigma(0)*(1+dlumi_lumi) + sigma(1)+sigmabar(1,:) -resum_expanded(1,:) 

    if (order >= 2) then
       matching_factor = matching_factor &
            &          + sigma(2) + sigmabar(2,:) - resum_expanded(2,:) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:))
    end if

    if (order >= 3) then
       matching_factor = matching_factor &
            &          + sigma(3) + sigmabar(3,:) - resum_expanded(3,:) &
            &          - (resum_expanded(2,:)/sigma(0)-dlumi_lumi_expanded(2,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigma(2)+sigmabar(2,:) - resum_expanded(2,:)) &
            &          + resum_expanded(1,:)/sigma(0)*(sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &                           * (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))
    end if

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor

  end function matching_a
  

  function matching_b(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    matching_factor = sigma(0)*(1+dlumi_lumi) +sigmabar(1,:) - resum_expanded(1,:)

    if (order >= 2) then
       matching_factor = matching_factor + sigma(1) + sigmabar(2,:) - resum_expanded(2,:)&
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &          * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:))
    end if

    if (order >= 3) then
       matching_factor = matching_factor &
            &          + sigma(2) + sigmabar(3,:) - resum_expanded(3,:) &
            &          - (resum_expanded(2,:)/sigma(0)-dlumi_lumi_expanded(2,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigma(2)+sigmabar(2,:) - resum_expanded(2,:)) &
            &          + resum_expanded(1,:)/sigma(0)*(sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &                           * (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))
    end if

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor

  end function matching_b


  function matching_c(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    if (order < 2) call wae_error("matching_c: scheme (c) not supported for fixed order < 2")

    if (order >= 2) then
       matching_factor = sigma(0)*(1+dlumi_lumi) + sigmabar(1,:) - resum_expanded(1,:) &
            &          + sigmabar(2,:) - resum_expanded(2,:) - sigma(1)/sigma(0)*sigmabar(1,:)&
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &          * (sigmabar(1,:) - resum_expanded(1,:))
    end if

    if (order >= 3) then
       matching_factor = matching_factor &
            &          + sigma(1) + sigmabar(3,:) - resum_expanded(3,:) &
            &          + sigma(1)/sigma(0)*sigmabar(1,:) - sigma(2)/sigma(0)*sigmabar(1,:) &
            &          - (resum_expanded(2,:)/sigma(0)-dlumi_lumi_expanded(2,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))*sigma(1) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigmabar(2,:) - resum_expanded(2,:)) &
            &          + resum_expanded(1,:)/sigma(0)*(sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &                           * (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))
    end if

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor

  end function matching_c


  ! the following routine implements matching c with the NLO K factor unexpanded
  ! it differs from matching c at N3LO
  function matching_c_unexp(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    if (order < 2) call wae_error("matching_c: scheme (c) not supported for fixed order < 2")

    if (order >= 2) then
       matching_factor = sigma(0)*(1+dlumi_lumi) + sigmabar(1,:) - resum_expanded(1,:) &
            &          + sigmabar(2,:) - resum_expanded(2,:) - sigma(1)/sigma(0)*sigmabar(1,:)&
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &          * (sigmabar(1,:) - resum_expanded(1,:))
    end if

    if (order >= 3) then
       matching_factor = matching_factor &
            &          + sigma(1) + sigmabar(3,:) - resum_expanded(3,:) &
            &          + sigma(1)/sigma(0)*sigmabar(1,:) - sigma(2)/(sigma(0)+sigma(1))*sigmabar(1,:) &
            &          - (resum_expanded(2,:)/sigma(0)-dlumi_lumi_expanded(2,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))*sigma(1) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigmabar(2,:) - resum_expanded(2,:)) &
            &          + resum_expanded(1,:)/sigma(0)*(sigma(1)+sigmabar(1,:) - resum_expanded(1,:)) &
            &                           * (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))
    end if

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor

  end function matching_c_unexp



  function matching_d(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    if (order < 3) call wae_error("matching_d: scheme (d) not supported for fixed order < 3")

    if (order >= 3) then
       matching_factor = sigma(0)*(1+dlumi_lumi) + sigmabar(1,:) - resum_expanded(1,:) &
            &          + sigmabar(2,:) - resum_expanded(2,:) - sigma(1)/sigma(0)*sigmabar(1,:)&
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &          * (sigmabar(1,:) - resum_expanded(1,:)) &
            &          + sigmabar(3,:) - resum_expanded(3,:) - sigma(1)/sigma(0)*sigmabar(2,:) &
            &          + (sigma(1)**2/sigma(0)**2 - sigma(2)/sigma(0))*sigmabar(1,:) &
            &          - (resum_expanded(2,:)/sigma(0)-dlumi_lumi_expanded(2,:))&
            &                           * (sigmabar(1,:) - resum_expanded(1,:)) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigmabar(2,:) - resum_expanded(2,:) - sigma(1)/sigma(0)*sigmabar(1,:)) &
            &          + resum_expanded(1,:)/sigma(0)*(sigmabar(1,:) - resum_expanded(1,:)) &
            &                           * (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))

       !matching_factor = sigma(0)*(one+dlumi_lumi) + sigmabar(1,:) - resum_expanded(1,:) &
       !     & + sigmabar(2,:) - resum_expanded(2,:) - sigma(1)/sigma(0)*sigmabar(1,:) &
       !     & + sigmabar(3,:) - resum_expanded(3,:) &
       !     & + (sigma(1)**2/sigma(0)**2-sigma(2)/sigma(0))*sigmabar(1,:) &
       !     & - sigma(1)/sigma(0)*sigmabar(2,:) &
       !     & + (dlumi_lumi_expanded(1,:) - resum_expanded(1,:)/sigma(0) &
       !     & + dlumi_lumi_expanded(2,:) - resum_expanded(2,:)/sigma(0))*(sigmabar(1,:) - resum_expanded(1,:)) &
       !     & + (dlumi_lumi_expanded(1,:) - resum_expanded(1,:)/sigma(0))*(sigmabar(2,:) &
       !     & - sigma(1)/sigma(0)*sigmabar(1,:) - resum_expanded(2,:)) &
       !     & - resum_expanded(1,:)/sigma(0)*(dlumi_lumi_expanded(1,:) &
       !     & - resum_expanded(1,:)/sigma(0))*(sigmabar(1,:) - resum_expanded(1,:))
    end if

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor

  end function matching_d


  ! THE FOLLOWING FUNCTION IS NOT TESTED AND NOT INTENDED FOR ANYTHING OTHER
  ! THAN INTERNAL USE BY THE AUTHORS
  function matching_logR(sigma, pt, resummation, resum_expanded, sigmabar) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    real(dp) :: Z(size(pt))

    ! the brackets look wrong for the sigma(0) and sigma(1) contributions 
    res = (resummation/sigma(0)+sigma(1))*exp((sigma(1)+sigmabar(1,:) -resum_expanded(1,:))/sigma(0)) &
         &                           *exp((sigmabar(2,:)-resum_expanded(2,:))/sigma(0) - &
         &((sigma(1)+sigmabar(1,:))**2 -resum_expanded(1,:)**2)/(two*sigma(0)**2))


  end function matching_logR


  function matching_R(sigma, pt, resummation, resum_expanded, sigmabar) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    real(dp) :: Z(size(pt)), cutoff
    integer  :: i, order

    order = size(resum_expanded,dim=1)
    ! pt cutoff at which we switch off the finite remainder (it can be set to zero)
    !cutoff = 2.0_dp
    cutoff = zero

    ! in this formulation, the whole two-loop constant term 
    ! is NOT suppressed by the Sudakov form factor
    res = resummation
     
    do i = 1, size(pt)
       if (pt(i) .ge. cutoff) then
          res(i) = res(i) + sigma(1)+sigmabar(1,i)-resum_expanded(1,i)
          if (order .ge. 2) res(i) = res(i) + sigma(2)+sigmabar(2,i)&
               & -resum_expanded(2,i)
          if (order .eq. 3) res(i) = res(i) + sigma(3)+sigmabar(3,i)&
               & -resum_expanded(3,i)
       end if
    end do
    

  end function matching_R


  ! THE FOLLOWING THREE MATCHING SCHEMES ARE DEFINED TO INVESTIGATE
  ! THE EFFECT OF HEAVY QUARK MASSES AT HIGHER ORDERS. THEY ARE NOT
  ! INTENDED FOR ANYTHING OTHER THAN INTERNAL USE BY THE AUTHORS
  function matching_moda(sigma, pt, cs, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    type(process_and_parameters), intent(in) :: cs
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))
    real(dp)             :: lNLL(size(pt)), L_tilde(size(pt)), sigma_native(0:1)

    order = size(resum_expanded,dim=1)

    L_tilde = Ltilde(pt/cs%Q, cs%p)
    sigma_native = cross_sections(cs)
    lNLL  = sigma(0)/sigma_native(0)*lumi_NLL(exp(-L_tilde)*cs%muF, cs)
    matching_factor = lNLL*(1+dlumi_lumi) + sigma(1)+sigmabar(1,:) -resum_expanded(1,:) 

    if (order >= 2) then
       matching_factor = matching_factor &
            &          + sigma(2) + sigmabar(2,:) - resum_expanded(2,:) &
            &          - (resum_expanded(1,:)/sigma(0)-dlumi_lumi_expanded(1,:))&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:))
    end if

    if (order >= 3) call wae_error("matching_moda: moda-scheme not supported for fixed order >= 3")

    res =  (resummation/lNLL/(1+dlumi_lumi))*matching_factor
  end function matching_moda


  function matching_modRa(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))

    order = size(resum_expanded,dim=1)

    matching_factor = sigma(0)*(1+dlumi_lumi)

    if (order >= 2) then
       matching_factor = matching_factor &
            &          + sigma(2) + sigmabar(2,:) - resum_expanded(2,:)
    end if

    if (order >= 3) call wae_error("matching_modR: modR-scheme not supported for fixed order >= 3")

    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor &
         & + sigma(1)+sigmabar(1,:) -resum_expanded(1,:)

  end function matching_modRa


  function matching_loga(sigma, pt, resummation, resum_expanded, sigmabar,dlumi_lumi,dlumi_lumi_expanded) result(res)
    real(dp), intent(in) :: sigma(0:), pt(:), resummation(:), resum_expanded(1:,:), sigmabar(1:,:)
    real(dp), intent(in) :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
    real(dp)             :: res(size(pt))
    !----------------------------------------------
    integer              :: order
    real(dp)             :: matching_factor(size(resummation))
    real(dp)             :: lumil10, norm

    order = size(resum_expanded,dim=1)

    matching_factor = sigma(0)*(1+dlumi_lumi)

    if (order >= 2) then
       matching_factor = matching_factor &
            &          + sigma(2) + sigmabar(2,:) - resum_expanded(2,:) &
            &          - resum_expanded(1,:)/sigma(0)&
            &                           * (sigma(1)+sigmabar(1,:) - resum_expanded(1,:))
       matching_factor = matching_factor &
            & *exp(-half*(sigma(1)+sigmabar(1,:) -resum_expanded(1,:))**2/sigma(0)**2)
      
    end if
    
    res =  (resummation/sigma(0)/(1+dlumi_lumi))*matching_factor
    res = res*exp((sigma(1)+sigmabar(1,:) -resum_expanded(1,:))/sigma(0))
    
    if (order > 2) then
       call wae_error("matching_loga: loga matching scheme is supported only for fixed_order = 1, 2")
    end if

    ! correct normalization for loga matching scheme (order 2)
    !lumil10 = dlumi_lumi_expanded(1,1)
    !norm = (sigma(0)*(one+lumil10)+sigma(2) &
    !     & -lumil10*(sigma(1)-sigma(0)*lumil10)) &
    !     & * exp(sigma(1)/sigma(0)-lumil10) &
    !     & * exp(-half*(sigma(1)/sigma(0)-lumil10)**2)
    !write(*,*) "norm = ", norm
   
  end function matching_loga


  !----------------------------------------------------------------------
  ! Function that returns the result in a given fixed-order scheme;
  ! the order to be used is deduced from the number of columns in the
  ! sigmabar array.
  !
  ! Input parameters:
  !   - scheme = string indicating the name of the scheme
  !   - sigma(0:fixed_order) = terms in the total cross section from LO(0) up to fixed_order
  !   - sigmabar(1:fixed_order, :) = 1-jet cross section, integrated from infinity down to pt (2nd index)
  !
  function fixed_order_schemes(scheme, sigma, sigmabar) result(res)
    character(len=*), intent(in) :: scheme
    real(dp), intent(in) :: sigma(0:), sigmabar(1:,:)
    real(dp)             :: res(size(sigmabar,dim=2)), norm
    integer              :: fixed_order 

    fixed_order = size(sigmabar,dim=1)
    if (size(sigma)-1 < fixed_order) call wae_error("fixed_order_schemes", &
         & "size of sigma(:) array is smaller than fixed_order deduced from sigmabar")

    select case(trim(scheme))
    case("a")
       norm = sum(sigma(:fixed_order))
       res  = one + sum(sigmabar,dim=1)/norm
    case("b")
       norm = sum(sigma(0:fixed_order-1))
       res  = one + sum(sigmabar,dim=1)/norm
    case("c")
       if (fixed_order < 2) call wae_error("fixed_order_schemes: scheme (c) unsupported for fixed_order < 2")
       norm = sum(sigma(0:fixed_order-2))
       res = one + sum(sigmabar,dim=1)/norm &
            &      - sigma(fixed_order-1)/norm/sigma(0)*sigmabar(1,:)
    case("cunx")
       ! Second version of the c scheme at N3LO
       ! with NLO K factor unexpanded
       norm = sum(sigma(0:fixed_order-2))
       res = one + sum(sigmabar,dim=1)/norm &
            &      - sigma(fixed_order-1)/norm**2*sigmabar(1,:)       
    case("d")
       ! check out the following expression
       if (fixed_order < 3) call wae_error("fixed_order_schemes: scheme (d) unsupported for fixed_order < 3")
       norm = sum(sigma(0:fixed_order-3))
       res = one + sum(sigmabar,dim=1)/norm &
            &      - sigma(fixed_order-2)/norm**2*(sigmabar(1,:)+sigmabar(2,:)) &
            &      + (sigma(1)*sigma(fixed_order-2)-sigma(0)*sigma(fixed_order-1))*sigmabar(1,:)/norm**3
    case("moda")
       norm = sum(sigma(:fixed_order))
       res  = one + sum(sigmabar,dim=1)/norm
    case("modRa")
       norm = sum(sigma(:fixed_order))
       res  = one + sum(sigmabar,dim=1)/norm
    case("R")
       norm = sum(sigma(:fixed_order))
       res  = one + sum(sigmabar,dim=1)/norm
    case("loga")
       norm = sum(sigma(:fixed_order))
       res  = one + sum(sigmabar,dim=1)/norm
    case default
       call wae_error("fixed_order_schemes: unrecognised fixed order scheme "//scheme)
    end select
  end function fixed_order_schemes
  
  
end module matching
