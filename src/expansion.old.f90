module expansion
  use types; use consts_dp
  use rad_tools; use emsn_tools
  use pdfs_tools
  use hoppet_v1
  implicit none

  private
  public :: expanded_sigma

contains

  !======================================================================
  function expanded_sigma(pt, cs, order,dlumi_lumi_expanded) result(res)
    real(dp),                  intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order !(0=LL,1=NLL,2=NNLL)
    real(dp), optional, intent(out)       :: dlumi_lumi_expanded(size(pt))
    ! first index of res denotes the order in alphas 
    real(dp)                              :: res(2,size(pt))
    !-------------------------------------
    real(dp) :: L(size(pt))
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: P0_x_pdf1(0:grid%ny,-6:7), P0_x_pdf2(0:grid%ny,-6:7)
    real(dp) :: lumi_LO_LO, lumi_P0_LO, lumi_LO_P0, lumi_LO_P02, lumi_P02_LO
    real(dp) :: lumi_LO_P1, lumi_P1_LO, lumi_P0_P0
    real(dp) :: lumi_C_plus_P_LO, lumi_LO_C_plus_P, lumi_C_plus_P_LO_P0, lumi_P0_C_plus_P_LO
    real(dp) :: lumi_C_plus_P_P0_LO, lumi_LO_C_plus_P_P0, lumi_C_plus_P_C_plus_P
    real(dp) :: G12(size(pt)), G11(size(pt))
    real(dp) :: mtwoL(size(L))
    real(dp) :: Log_muR_M, Log_M_Q, x, Log_muF_muR
    type(split_mat) :: C_plus_P

    if (present(dlumi_lumi_expanded)) dlumi_lumi_expanded = zero

    L = Ltilde(pt/cs%Q, cs%p)

    ! prepare the material needed for the PDFs
    call get_pdfs(cs%muF, cs%collider, pdf1, pdf2)

    lumi_LO_LO = lumi_at_x(pdf1, pdf2, cs)

    if (order >= order_NLL) then
       P0_x_pdf1   = dglap_h%P_LO .conv. pdf1
       P0_x_pdf2   = dglap_h%P_LO .conv. pdf2

       lumi_P0_LO  = lumi_at_x(P0_x_pdf1, pdf2, cs)/twopi
       lumi_LO_P0  = lumi_at_x(pdf1, P0_x_pdf2, cs)/twopi

       lumi_P02_LO = lumi_at_x(dglap_h%P_LO .conv. P0_x_pdf1, pdf2, cs)/twopi**2
       lumi_LO_P02 = lumi_at_x(pdf1, dglap_h%P_LO .conv. P0_x_pdf2, cs)/twopi**2

       lumi_P1_LO  = lumi_at_x(dglap_h%P_NLO .conv. pdf1, pdf2, cs)/twopi**2
       lumi_LO_P1  = lumi_at_x(pdf1, dglap_h%P_NLO .conv. pdf2, cs)/twopi**2

       lumi_P0_P0  = lumi_at_x(P0_x_pdf1, P0_x_pdf2, cs)/twopi**2

    end if

    if (order >= order_NNLL) then

       call InitSplitMat(C_plus_P, C1_matrix)
       call AddWithCoeff(C_plus_P, dglap_h%P_LO, -cs%ln_muF2_M2+cs%ln_Q2_M2)

       lumi_C_plus_P_LO    = lumi_at_x(C_plus_P .conv. pdf1, pdf2, cs)/twopi
       lumi_LO_C_plus_P    = lumi_at_x(pdf1, C_plus_P .conv. pdf2, cs)/twopi

       lumi_C_plus_P_LO_P0 = lumi_at_x(C_plus_P .conv. pdf1, dglap_h%P_LO .conv. pdf2, cs)/twopi**2
       lumi_P0_C_plus_P_LO = lumi_at_x(dglap_h%P_LO .conv. pdf1, C_plus_P .conv. pdf2, cs)/twopi**2

       lumi_C_plus_P_P0_LO = lumi_at_x(C_plus_P .conv. (dglap_h%P_LO .conv. pdf1), pdf2, cs)/twopi**2
       lumi_LO_C_plus_P_P0 = lumi_at_x(pdf1, C_plus_P .conv. (dglap_h%P_LO .conv. pdf2), cs)/twopi**2

       lumi_C_plus_P_C_plus_P = lumi_at_x(C_plus_P .conv. pdf1, C_plus_P .conv. pdf2, cs)/twopi**2

    end if

    res = zero

    !-- LL terms
    G12 = -2*CC*L**2/Pi
    res(1,:) = G12 * lumi_LO_LO
    res(2,:) = ((2*CC**2*L**4)/Pi**2 -8*beta0*CC*L**3/(3._dp*Pi)) * lumi_LO_LO 

    !-- NLL terms
    if (order >= order_NLL) then
       G11 = (-4*(BB*CC*L + CC*L*(-half*cs%ln_Q2_M2) ))/Pi
       res(1,:) = res(1,:) + lumi_LO_LO * G11

       Log_muR_M   =  half*cs%ln_muR2_M2
       Log_M_Q     = -half*cs%ln_Q2_M2
       Log_muF_muR =  half*(cs%ln_muF2_M2 - cs%ln_muR2_M2)

       res(2,:) = res(2,:) + lumi_LO_LO * (8*BB**2*CC**2*L**2 - CC*cmw_K*L**2 + 8*BB*CC**2*L**3 &
            &              - 4*BB*beta0*CC*L**2*Pi - 4*beta0*CC*L**2*Pi*Log_muR_M  &
            &              + 16*BB*CC**2*L**2*Log_M_Q + 8*CC**2*L**3*Log_M_Q &
            &              - 8*beta0*CC*L**2*Pi*Log_M_Q + 8*CC**2*L**2*Log_M_Q**2)/Pi**2

       !----
       mtwoL = -two*L ! the evolution distance in ln Q^2
       res(1,:) = res(1,:) + mtwoL * (lumi_P0_LO + lumi_LO_P0)
       res(2,:) = res(2,:) + (G12+G11) * (mtwoL) * (lumi_P0_LO + lumi_LO_P0) &
            &              + mtwoL**2 * (half*(lumi_P02_LO + lumi_LO_P02) &
            &              + lumi_P0_P0 - half*beta0*(lumi_P0_LO + lumi_LO_P0)) & 
            &              + mtwoL * (lumi_P1_LO + lumi_LO_P1 &
            &              + beta0 * (cs%ln_muR2_M2 - cs%ln_muF2_M2)*(lumi_P0_LO + lumi_LO_P0)) 
    end if


    ! NNLL terms
    if (order >= order_NNLL) then

       res(1,:) = res(1,:) + lumi_LO_LO*H(1)/twopi
       res(2,:) = res(2,:) - lumi_LO_LO*CC*H(1)/pisq*L**2 &
            &              + lumi_LO_LO*(-B(2)/twopi/pi - 2._dp*BB*CC*H(1)/pisq  &
            &              - 2._dp*CC*H(1)/pisq*Log_M_Q - 8._dp*BB*CC*beta0/pi*(Log_muR_M+Log_M_Q) &
            &              - 8._dp*CC*beta0/pi*Log_M_Q*(Log_muR_M+Log_M_Q) - 2._dp*CC*cmw_K/pisq*Log_M_Q &
            &              + 4._dp*CC/pisq*non_incl(cs%jet_radius,'all'))*L &
            &              - H(1)/pi*(lumi_P0_LO + lumi_LO_P0)*L 

       !----    
       res(1,:) = res(1,:) + lumi_LO_C_plus_P + lumi_C_plus_P_LO 
       res(2,:) = res(2,:) - 2._dp*CC/pi*(lumi_LO_C_plus_P + lumi_C_plus_P_LO)*L**2    &
            &              - 4._dp*BB*CC/pi*(lumi_LO_C_plus_P + lumi_C_plus_P_LO)*L      &
            &              - 4._dp*CC/pi*Log_M_Q*(lumi_LO_C_plus_P + lumi_C_plus_P_LO)*L &
            &              + 2._dp*beta0*(lumi_LO_C_plus_P + lumi_C_plus_P_LO)*L         &
            &              - 2._dp*(lumi_C_plus_P_LO_P0 + lumi_P0_C_plus_P_LO)*L         &
            &              - 2._dp*(lumi_LO_C_plus_P_P0 + lumi_C_plus_P_P0_LO)*L         

       ! Spurious alphas^2 constant term ---
       if (cs%include_c1_squared) res(2,:) = res(2,:) &
            &              + H(1)/twopi*(lumi_C_plus_P_LO + lumi_LO_C_plus_P) &
            &              + lumi_C_plus_P_C_plus_P                                 !&
       ! The following terms restore Q-invariance of the partonic cross section at alphas^2
       ! It MUST be used for debugging, NOT included in matching 
       !&              + lumi_LO_LO*non_incl(cs%jet_radius,'all')*A(1)*2._dp*(Log_muR_M+Log_M_Q)/pisq &
       !&              + lumi_LO_LO*H(2)/twopi**2

       if (present(dlumi_lumi_expanded)) then
          dlumi_lumi_expanded = H(1)/twopi + (lumi_C_plus_P_LO + lumi_LO_C_plus_P)/lumi_LO_LO
          dlumi_lumi_expanded = dlumi_lumi_expanded * cs%alphas_muR
       end if
    end if

    res(1,:) = res(1,:) * cs%alphas_muR
    res(2,:) = res(2,:) * cs%alphas_muR**2

  end function expanded_sigma


end module expansion
