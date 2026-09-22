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
    case("a","anew")
       ! 'a' is retained only as the internal polynomial-algebra identifier.
       res = matching_a(sigma, pt(:), resummation, resum_expanded, &
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
  

  function fixed_order_schemes(scheme, sigma, sigmabar) result(res)
    character(len=*), intent(in) :: scheme
    real(dp), intent(in) :: sigma(0:), sigmabar(1:,:)
    real(dp)             :: res(size(sigmabar,dim=2)), norm
    integer              :: fixed_order 

    fixed_order = size(sigmabar,dim=1)
    if (size(sigma)-1 < fixed_order) call wae_error("fixed_order_schemes", &
         & "size of sigma(:) array is smaller than fixed_order deduced from sigmabar")

    select case(trim(scheme))
    case("a","anew")
       norm = sum(sigma(:fixed_order))
       res = one + sum(sigmabar,dim=1)/norm
    case default
       call wae_error("fixed_order_schemes: unrecognised fixed order scheme "//scheme)
    end select
  end function fixed_order_schemes
  
  
end module matching
