module expansion
  use types; use consts_dp
  use rad_tools; use emsn_tools
  use pdfs_tools; use pdf_expansions
  use hoppet_v1 
  use ew_parameters; use mass_corr
  
  implicit none

  private
  public :: expanded_sigma

contains

  !======================================================================
  function expanded_sigma(pt, cs, order,matching_scheme,dlumi_lumi_expanded) result(res)
    real(dp),                  intent(in) :: pt(:)
    type(process_and_parameters), intent(in) :: cs
    integer,                   intent(in) :: order !(0=LL,1=NLL,2=NNLL)
    real(dp), optional, intent(out)       :: dlumi_lumi_expanded(2,size(pt))
    character(len=*), intent(in) :: matching_scheme
    ! first index of res denotes the order in alphas 
    real(dp)                              :: res(3,size(pt))
    !-------------------------------------
    real(dp) :: L(size(pt))
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: P0_x_pdf1(0:grid%ny,-6:7), P0_x_pdf2(0:grid%ny,-6:7)
    real(dp) :: G12(size(pt)), G11(size(pt))
    real(dp) :: mtwoL(size(L))
    real(dp) :: Log_muR_M, Log_M_Q, x, Log_muF_muR
    real(dp) :: h11, h12, h21, h22, h23, h24
    real(dp) :: h36, h35, h34, h33, h32, h31

    type(PDFExpansion)       :: pdf1e, pdf2e
    type(LumiExpansion)      :: lumie
    type(LumiExpansion_at_x) :: lumi
    real(dp)                 :: u 

    if (present(dlumi_lumi_expanded)) dlumi_lumi_expanded = zero

    L = Ltilde(pt/cs%Q, cs%p)
    
    ! prepare the material needed for the PDFs
    call get_pdfs(cs%muF, cs%collider, pdf1, pdf2)
!------------------------------------------------------------
!------------------------------------------------------------
    call AllocPDFExpansion(pdf1e)
    call AllocPDFExpansion(pdf2e)
    call AllocLumiExpansion(lumie)  
    !-- calculate expansion ---------------------------

    u = cs%ln_muF2_M2 - cs%ln_muR2_M2

    pdf1e%p = pdf1
    call set_pdf_expansion(pdf1e,u, cs)
    pdf2e%p = pdf2
    call set_pdf_expansion(pdf2e,u, cs)

!--- the following call returns luminosities factors lij, where lij is 
!--- the coefficient of (as2pi(muR))^i * t^j with  t = 2*log(muF1/muF0)=-2*L    
    if (order <= order_NLL) then
       call set_lumi_expansion(pdf1e,pdf2e,lumie, .false., .false., .false., cs) 
    endif
    if (order > order_NLL) then
       if (cs%include_c1_squared) then 
          call set_lumi_expansion(pdf1e,pdf2e,lumie, .true., .true., .true., cs)
       else
          call set_lumi_expansion(pdf1e,pdf2e,lumie, .true., .true., .false., cs)
       endif
    endif
 !--- evaluate luminosity at x=M^2/rts^2   
    call set_lumi_expansion_at_x(lumi,lumie,cs)
    

    res = zero

    !-- LL terms
    h12 = -two**2*CC
    h11 = zero
    h21 = zero
    h22 = zero
    h23 = -32._dp*CC*Pi*beta0/(3._dp)
    h24 = two**3*CC**2
    h36 = -32._dp/3._dp*CC**3
    h35 = 128._dp/3._dp*CC**2*beta0*Pi
    h34 = -32._dp*CC*beta0**2*Pi**2
    h33 = zero
    h32 = zero
    ! - Added by FD -
    h31 = zero

    !-- NLL terms
    if (order >= order_NLL) then

       h11 = h11 - 8._dp*BB*CC+4._dp*CC*cs%ln_Q2_M2
       h22 = h22 + 32*BB**2*CC**2-16*BB*CC*beta0*Pi-32*BB*CC**2*cs%ln_Q2_M2 &
            & + 16*CC*beta0*Pi*cs%ln_Q2_M2+8*CC**2*cs%ln_Q2_M2**2-8*CC*beta0*Pi*cs%ln_muR2_M2 &
            & - 4*CC*cmw_K
       h23 = h23 + 32*BB*CC**2-16*CC**2*cs%ln_Q2_M2 
       h35 = h35 - 64._dp*BB*CC**3+32._dp*CC**3*cs%ln_Q2_M2
       h34 = h34 -128*BB**2*CC**3 + (448*BB*beta0*CC**2*Pi)/3._dp + 128*BB*CC**3*cs%ln_Q2_M2 &
            & - (224*beta0*CC**2*Pi*cs%ln_Q2_M2)/3._dp - 32*CC**3*cs%ln_Q2_M2**2 & 
            & + 32*beta0*CC**2*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2) + 16*CC**2*cmw_K

       h33 = h33 + (-256*BB**3*CC**3)/3._dp + 128*BB**2*beta0*CC**2*Pi - (128*BB*beta0**2*CC*Pi**2)/3._dp & 
            & - (64*beta1*CC*Pi**2)/3._dp + 128*BB**2*CC**3*cs%ln_Q2_M2 &
            & - 128*BB*beta0*CC**2*Pi*cs%ln_Q2_M2 + (64*beta0**2*CC*Pi**2*cs%ln_Q2_M2)/3._dp &
            & - 64*BB*CC**3*cs%ln_Q2_M2**2 + 32*beta0*CC**2*Pi*cs%ln_Q2_M2**2 &
            & + (32*CC**3*cs%ln_Q2_M2**3)/3._dp +  64*BB*beta0*CC**2*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & - (128*beta0**2*CC*Pi**2*(cs%ln_muR2_M2-cs%ln_Q2_M2))/3._dp &
            & - 32*beta0*CC**2*Pi*cs%ln_Q2_M2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & + 32*BB*CC**2*cmw_K - (64*beta0*CC*Pi*cmw_K)/3._dp &
            & - 16*CC**2*cs%ln_Q2_M2*cmw_K
    end if


    ! NNLL terms
    if (order > order_NLL) then

       h21 = h21 - two*B(2)-16*BB*CC*beta0*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & +8*CC*beta0*Pi*cs%ln_Q2_M2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & +4._dp*CC*cmw_K*cs%ln_Q2_M2

       if (cs%observable == 'ptj') then 
          h21 = h21 + 16*CC*non_incl(cs%jet_radius,'all')
          ! - Added by FD -
          if (cs%small_r) then
             h31 = h31 + 16*CC*non_incl_aslnR_sq(cs%jet_radius,cs%small_r_R0)
          end if
       else
          h21 = h21 - 16*CC**2*zeta3
       endif

       h33 = h33 + 8*B(2)*CC + 64*BB*beta0*CC**2*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & - 32*beta0*CC**2*Pi*cs%ln_Q2_M2*(cs%ln_muR2_M2-cs%ln_Q2_M2)  &
            & - 16*CC**2*cs%ln_Q2_M2*cmw_K - 64*CC**2*non_incl(cs%jet_radius,'all')

       h32 = h32 - 2*A(3) + 16*B(2)*BB*CC - 8*B(2)*beta0*Pi - 32*BB*beta1*CC*Pi**2 &
            & - 8*B(2)*CC*cs%ln_Q2_M2 + 16*beta1*CC*Pi**2*cs%ln_Q2_M2 &
            & + 128*BB**2*beta0*CC**2*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & - 64*BB*beta0**2*CC*Pi**2*(cs%ln_muR2_M2-cs%ln_Q2_M2)  &
            & - 16*beta1*CC*Pi**2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & - 128*BB*beta0*CC**2*Pi*cs%ln_Q2_M2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & + 32*beta0**2*CC*Pi**2*cs%ln_Q2_M2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & + 32*beta0*CC**2*Pi*cs%ln_Q2_M2**2*(cs%ln_muR2_M2-cs%ln_Q2_M2) &
            & - 16*beta0**2*CC*Pi**2*(cs%ln_muR2_M2-cs%ln_Q2_M2)**2 &
            & - 32*BB*CC**2*cs%ln_Q2_M2*cmw_K + 16*beta0*CC*Pi*cs%ln_Q2_M2*cmw_K &
            & + 16*CC**2*cs%ln_Q2_M2**2*cmw_K - 16*beta0*CC*Pi*(cs%ln_muR2_M2-cs%ln_Q2_M2)*cmw_K &
            & - 128*BB*CC**2*non_incl(cs%jet_radius,'all') &
            & + 128*beta0*CC*Pi*non_incl(cs%jet_radius,'all') &
            & + 64*CC**2*cs%ln_Q2_M2*non_incl(cs%jet_radius,'all')

    end if

    !---------------------------------------------------------------------------------

    ! Build the expansion
    !-- LL terms
    res(1,:) = lumi%l00*h12*L**2
    res(2,:) = lumi%l00*(h24*L**4 + h23*L**3)
    res(3,:) = lumi%l00*(h36*L**6 + h35*L**5 + h34*L**4)
    !-- NLL and NNLL terms 
    if (order >= order_NLL) then
       res(1,:) = res(1,:) + (lumi%l00*h11 - two*lumi%l11)*L + lumi%l10
       if ((cs%loop_mass .ne. 'none') .and. (order > order_NLL)) then
          ! lumie%l10_mass-effects = lumie%l10 + (Hm(1)-H(1))*lumie%l00
          res(1,:) = res(1,:) + (Hm(1)-H(1))*lumi%l00
       endif
       res(2,:) = res(2,:) + (-two*lumi%l11*h12)*L**3 &
            & + (lumi%l00*h22 + lumi%l10*h12 - two*lumi%l11*h11 + four*lumi%l22)*L**2 &
            & + (lumi%l00*h21 + lumi%l10*h11 - two*lumi%l21)*L  + lumi%l20

! Internal checks only: the following (N3LL) piece appears only if the full 
! NNLL rp is used in the F[rp] function
!       if (cs%observable == 'ptB' .and. order == order_NNLL) then 
!          ! add subleading piece due to the hard-collinear terms in Rad_pNNLL
!          res(2,:) = res(2,:) - two*zeta3*A(1)*(B(1)+A(1)*(-cs%ln_Q2_M2))*lumi%l00
!       endif

       res(3,:) = res(3,:) + (-two*lumi%l11*h24)*L**5 &
            & + (lumi%l10*h24 - two*lumi%l11*h23 + four*lumi%l22*h12)*L**4 &
            & + (lumi%l00*h33 + lumi%l10*h23 - two*lumi%l11*h22 - two*lumi%l21*h12 &
            & + four*lumi%l22*h11 - two*four*lumi%l33)*L**3 &
            & + (lumi%l00*h32 + lumi%l10*h22 - two*lumi%l11*h21 + lumi%l20*h12 &
            & - two*lumi%l21*h11 + four*lumi%l32)*L**2 &
            & + (lumi%l00*h31 + lumi%l10*h21 + lumi%l20*h11 - two*lumi%l31)*L + lumi%l30
    end if
    
    ! rescale as^2 and as^3 terms according to hnnlo-v2.0 when heavy-quark
    ! mass effects are included (uncomment to use the EFT expansion beyond NLO)
    !if (cs%loop_mass .ne. 'none') then
    !    write(*,*) "Warning: mass effects are only included up to NLO"
    !   res(2,:) = res(2,:) * norm_mass(cs%M, mt, mb, 't')/norm_mass(cs%M, mt, mb, cs%loop_mass)
    !   res(3,:) = res(3,:) * norm_mass(cs%M, mt, mb, 't')/norm_mass(cs%M, mt, mb, cs%loop_mass)
    !endif
    

    if (present(dlumi_lumi_expanded)) then
       dlumi_lumi_expanded(1,:) = lumi%l10/lumi%l00 * (cs%alphas_muR/twopi)
       if (matching_scheme == 'moda') then
          dlumi_lumi_expanded(1,:) = dlumi_lumi_expanded(1,:) &
               & - two*L*lumi%l11/lumi%l00 * (cs%alphas_muR/twopi)
       end if
       if (cs%loop_mass .ne. 'none'.and.(order > order_NLL)) then
          ! lumie%l10_mass-effects = lumie%l10 + (Hm(1)-H(1))*lumie%l00
          dlumi_lumi_expanded(1,:) = dlumi_lumi_expanded(1,:) + (Hm(1)-H(1)) * (cs%alphas_muR/twopi)
       endif
       ! mass corrections to dlumi_lumi_expanded(2,:) are neglected since it
       ! only contributes to the N3LO matching!
       dlumi_lumi_expanded(2,:) = (-two*L*(lumi%dl21 - lumi%l10*lumi%l11/lumi%l00) &
            & + lumi%l20)/lumi%l00 * (cs%alphas_muR/twopi)**2
    end if


    res(1,:) = res(1,:) * (cs%alphas_muR/twopi)
    res(2,:) = res(2,:) * (cs%alphas_muR/twopi)**2
    res(3,:) = res(3,:) * (cs%alphas_muR/twopi)**3

    call Delete(pdf1e)
    call Delete(pdf2e)
    call Delete(lumie)

  end function expanded_sigma


end module expansion
