!========================================================
!--------------------------------------------------------
! Includes a module providing the parton luminosities
! and their convolutions with splitting and coefficient
! functions @ LO, NLO, NNLO
! -------------------------------------------------------
!========================================================
module pdfs_tools
  use types; use consts_dp
  !! if using LHAPDF, rename a couple of hoppet functions which
  !! would otherwise conflict with LHAPDF 
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF
  use ew_parameters
  use rad_tools

  implicit none

  type(dglap_holder), save, public :: dglap_h ! splitting function holder
  type(pdf_table) :: PDFs  ! parton densities
  type(grid_def), public :: grid
  type(split_mat), save, public :: C1_matrix

  private 
  public :: init_pdfs_from_LHAPDF
  public :: init_pdfs_hoppet_evolution
  public :: get_pdfs, lumi_at_x, lumi_all_x, lumi_LL, lumi_NLL, lumi_NNLL, dlumi_NLL 

  public :: cross_sections
contains

  !=======================================================================================
  ! set up Hoppet's x grid, splitting functions and PDF tabulation
  ! object (just the memory allocation -- it gets filled below).
  subroutine init_grid_and_dglap(Qmax)
    use coefficient_functions_VH 
    real(dp), intent(in), optional :: Qmax
    type(grid_def) ::  gdarray(4) ! grid
    integer :: pdf_fit_order, nf_lcl
    real(dp) :: dy, ymax, dlnlnQ
    ! build the PDF grid 
    dy = 0.10_dp
    dlnlnQ = dy/4.0_dp
    ymax = 15.0_dp
    pdf_fit_order = -5 
    call InitGridDef(gdarray(4),dy/27.0_dp, 0.2_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(3),dy/9.0_dp,  0.5_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(2),dy/3.0_dp,  2.0_dp, order=pdf_fit_order)
    call InitGridDef(gdarray(1),dy,         ymax  , order=pdf_fit_order)
    call InitGridDef(grid,gdarray(1:4),locked=.true.)


    nf_lcl = 5
    call qcd_SetNf(nf_lcl)
    call InitDglapHolder(grid, dglap_h,  factscheme = factscheme_MSbar, nloop = 3)
    call InitCoeffMatrix(grid, C1_matrix)

    if (present(Qmax)) then
       call AllocPdfTable(grid, PDFs, 1d0, Qmax, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    else
       call AllocPdfTable(grid, PDFs, 1d0, 2d4, dlnlnQ=dlnlnQ, freeze_at_Qmin = .true.)
    end if

  end subroutine init_grid_and_dglap

  !=======================================================================================
  !! initialise Hoppet's PDF table from LHAPDF
  subroutine init_pdfs_from_LHAPDF(pdf_name, pdf_set)
    character(len=*), intent(in) :: pdf_name
    integer,          intent(in) :: pdf_set
    !--------------------------------------------------
    real(dp) :: alphasMZ, alphasPDF
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    call init_grid_and_dglap()

    !------------------------------------------------------------
    ! set up LHAPDF
    call InitPDFsetByName(trim(pdf_name))
    call InitPDF(pdf_set)

    ! sort out the coupling: get it from LHAPDF, which will have been initialised in init_pdfs
    alphasMZ = alphasPDF(MZ)  
    call InitRunningCoupling(alfas=alphasMZ, Q=MZ, nloop=3, fixnf=nf_int)
    
    call FillPdfTable_LHAPDF(PDFs, evolvePDF)

  end subroutine init_pdfs_from_LHAPDF

  !=======================================================================================
  !! initialise Hoppet's PDF table by evolution from an LHAPDF input
  !! at some reference scale.
  subroutine init_pdfs_hoppet_evolution(pdf_name, pdf_set, ref_Q, alphas_Q, rts)
    character(len=*), intent(in) :: pdf_name
    integer,          intent(in) :: pdf_set
    real(dp),         intent(in) :: ref_Q, alphas_Q, rts
    !--------------------------------------------------
    real(dp), pointer :: pdf_Q(:,:)
    real(dp) :: Qmax
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    Qmax = max(2e4_dp, rts)
    call init_grid_and_dglap(rts)
    call AllocPDF(grid, pdf_Q)

    !------------------------------------------------------------
    ! set up LHAPDF
    if (pdf_name == "dummy") then
       pdf_Q = unpolarized_dummy_pdf(xValues(grid))
    else
       call InitPDFsetByName(trim(pdf_name))
       call InitPDF(pdf_set)
       ! extract the PDF at our reference Q value
       call InitPDF_LHAPDF(grid, pdf_Q, evolvePDF, ref_Q)
    end if

    ! set up a running coupling, set up the PDF tabulation and perform the evolution
    call InitRunningCoupling(alfas=alphas_Q, Q=ref_Q, nloop=3, fixnf=nf_int)
    call EvolvePdfTable(PDFs, ref_Q, pdf_Q, dglap_h, ash_global)
    !------------------------------------------------------------
  end subroutine init_pdfs_hoppet_evolution

  !======================================================================
  ! this works out the partonic luminosity for the LO of process proc, 
  ! returns as an array over values of x = M^2/rts^2; gives a luminosity
  ! based on the two input PDFs
  function lumi_all_x(pdf1, pdf2, cs) result(res)
    real(dp),                     intent(in) :: pdf1(:,-6:), pdf2(:,-6:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                     :: res(size(pdf1,1))
    !----------------------------
    integer i

    select case(trim(cs%proc))
    case("H")
       res = PartonLuminosity(grid, pdf1(:,iflv_g), pdf2(:,iflv_g))
    case('DY')
       res = zero
       do i = 1, nf_int
          res = res + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:, i), pdf2(:,-i))
          res = res + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:,-i), pdf2(:, i))
       end do
    end select
    res = res * cs%norm_fb ! Normalise in order to get the cross section in [fb]
  end function lumi_all_x

  !======================================================================
  ! this works out the partonic luminosity for the LO of process proc, 
  ! returns at a fixed value x = M^2/rts^2; gives a luminosity
  ! based on the two input PDFs
  function lumi_at_x(pdf1, pdf2, cs) result(res)
    real(dp),                     intent(in) :: pdf1(:,-6:), pdf2(:,-6:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                     :: res
    !--------------------------------------
    res = lumi_all_x(pdf1, pdf2, cs) .atx. (cs%M2_rts2 .with. grid)
  end function lumi_at_x

  !======================================================================
  ! returns the luminosity that enters the cross section at LL, i.e.
  ! evaluated at the default factorization scale for the process;
  function lumi_LL(cs) result(res)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                                 :: res
    !----------------------------------------------
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    integer  :: i

    call get_pdfs(cs%muF, cs%collider, pdf1, pdf2)
    res = lumi_at_x(pdf1, pdf2, cs)
  end function lumi_LL

  !======================================================================
  ! returns the luminosity that enters the cross section at NLL, which
  ! is to be evaluated separately for each of the input muF values
  ! (corresponding to different values of ptveto)
  !
  ! Eq.(25) of BMSZ
  function lumi_NLL(muF, cs) result(res)
    real(dp),                  intent(in) :: muF(:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                              :: res(size(muF))
    !----------------------------------------------
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    integer  :: i

    do i = 1, size(muF)
       call get_pdfs(muF(i), cs%collider, pdf1, pdf2)
       res(i) = lumi_at_x(pdf1, pdf2, cs)
    end do
  end function lumi_NLL


  !======================================================================
  ! returns the Eq. (42) of 1206.4998 
  function dlumi_NLL(muF, cs) result(res)
    real(dp),                  intent(in) :: muF(:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                              :: res(size(muF))
    !----------------------------------------------
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    integer  :: i

    do i = 1, size(muF)
       call get_pdfs(muF(i), cs%collider, pdf1, pdf2)
       res(i) = lumi_at_x((dglap_h%P_LO .conv. pdf1), pdf2, cs)
       res(i) = res(i)+ lumi_at_x(pdf1, (dglap_h%P_LO .conv. pdf2), cs)
    end do
    res = res * two * cs%as2pi
  end function dlumi_NLL


  !======================================================================
  ! returns the luminosity that enters the cross section at NNLL,
  ! including C1 terms. It is evaluated at each of the supplied muF
  ! values, and the user must supply to corresponding lambda values,
  ! which are of relevance for establishing the coupling that
  ! multiplies the C1 pieces.
  !
  ! Eq.(25)+Eq.(26) of BMSZ, but without the H term
  function lumi_NNLL(muF, lambda, cs) result(res)
    real(dp),                  intent(in) :: muF(:), lambda(:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                              :: res(size(muF))
    !----------------------------------------------
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: as2pi_mu
    type(split_mat) :: C_plus_P
    integer  :: i
    type(pdf_rep) :: pdfrep

    call InitSplitMat(C_plus_P, C1_matrix)
    call AddWithCoeff(C_plus_P, dglap_h%P_LO, -cs%ln_muF2_M2+cs%ln_Q2_M2)

    Do i = 1, size(muF)
       call get_pdfs(muF(i), cs%collider, pdf1, pdf2)
       ! get the running coupling at scale muR = muR_default * exp(-L)
       as2pi_mu = cs%as2pi / (1-2*lambda(i))

       if (cs%include_c1_squared) then
          ! Option 1: correct the PDFs and get an implicit as2pi_mu^2 term
          pdf1 = pdf1 + as2pi_mu * (C_plus_P .conv. pdf1)

          pdf2 = pdf2 + as2pi_mu * (C_plus_P .conv. pdf2)
          res(i) = lumi_at_x(pdf1, pdf2, cs)
       else 
          ! option 2: correct the luminosities and get no "spurious" as2pi_mu^2 term
          res(i) = lumi_at_x(pdf1, pdf2, cs)
          res(i) = res(i) + as2pi_mu*lumi_at_x(C_plus_P .conv. pdf1, pdf2, cs)
          res(i) = res(i) + as2pi_mu*lumi_at_x(pdf1, C_plus_P .conv. pdf2, cs)
       end if
    end do
  end function lumi_NNLL

  !======================================================================
  ! return the two incoming PDFs evaluated at scale muF, taking into
  ! account whether the collider is pp or ppbar
  subroutine get_pdfs(muF, collider, pdf1, pdf2)
    real(dp),         intent(in)  :: muF
    character(len=*), intent(in)  :: collider
    real(dp),         intent(out) :: pdf1(:,-6:), pdf2(:,-6:)
    !---------------------------------------------------------

    call EvalPdfTable_Q(PDFs, muF, pdf1)
    select case(trim(collider))
    case("pp")
       pdf2 = pdf1
    case("ppbar")
       pdf2(:,-6:6) = pdf1(:,6:-6:-1)
       pdf2(:,7)    = pdf1(:,7)          ! index 7 contains info on representation in flavour space
    case default
       call wae_error("unrecognized collider: "//collider)
    end select
  end subroutine get_pdfs


  !======================================================================
  !! The dummy PDF suggested by Vogt as the initial condition for the 
  !! unpolarized evolution (as used in hep-ph/0511119).
  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals),ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))
    !---------------------
    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv=5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half*N_ls

    pdf = zero
    ! clean method for labelling as PDF as being in the human representation
    ! (not actually needed after setting pdf=0
    call LabelPdfAsHuman(pdf)

    !-- remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1-xvals)**3
    dv = N_dv * xvals**0.8_dp * (1-xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1-xvals)**6
    ubar = dbar * (1-xvals)

    ! labels iflv_g, etc., come from the hoppet_v1 module, inherited
    ! from the main program
    pdf(:, iflv_g) = N_g * xvals**(-0.1_dp) * (1-xvals)**5
    pdf(:,-iflv_s) = 0.2_dp*(dbar + ubar)
    pdf(:, iflv_s) = pdf(:,-iflv_s)
    pdf(:, iflv_u) = uv + ubar
    pdf(:,-iflv_u) = ubar
    pdf(:, iflv_d) = dv + dbar
    pdf(:,-iflv_d) = dbar
  end function unpolarized_dummy_pdf


  !======================================================================
  ! returns an array containing the LO and NLO cross sections 
  function cross_sections(cs) result(res)
    use coefficient_functions_VH 
    type(process_and_parameters), intent(in) :: cs
    real(dp)                              :: res(1:2)
    !----------------------------------------------
    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: P_x_pdf1(0:grid%ny,-6:7), P_x_pdf2(0:grid%ny,-6:7)
    real(dp) :: lumi_gg(0:grid%ny), lumi_qg(0:grid%ny), lumi_qqbar(0:grid%ny)
    integer  :: i
    type(gdval) :: x

    call get_pdfs(cs%muF, cs%collider, pdf1, pdf2)
    x = cs%M2_rts2 .with. grid

    select case(cs%proc)
    case("H")
       call InitHiggsCoeffs(grid)
       lumi_gg = PartonLuminosity(grid, pdf1(:,iflv_g), pdf2(:,iflv_g))
       lumi_qg = PartonLuminosity(grid, pdf1(:,iflv_g), sum(pdf2(:,-6:-1),dim=2)&
            &                                          +sum(pdf2(:, 1:6 ),dim=2) )
       lumi_qg = lumi_qg + &
            &    PartonLuminosity(grid, sum(pdf1(:,-6:-1),dim=2)&
            &                          +sum(pdf1(:, 1:6 ),dim=2), pdf2(:,iflv_g)   )
       lumi_qqbar = 0
       do i = -6, 6
          if (i == 0) cycle
          lumi_qqbar = lumi_qqbar + PartonLuminosity(grid, pdf1(:,i), pdf2(:,-i))
       end do

       res(1) = lumi_gg .atx. x

       res(2) =   ((higgs_H1gg .conv. lumi_gg) .atx. x)&
            &   + ((higgs_H1qg .conv. lumi_qg) .atx. x)&
            &   + ((higgs_H1qqbar .conv. lumi_qqbar) .atx. x)

       if (cs%ln_muF2_M2 /= zero) then
          P_x_pdf1 = dglap_h%P_LO * pdf1
          P_x_pdf2 = dglap_h%P_LO * pdf2
          res(2) = res(2) - (  (PartonLuminosity(grid, P_x_pdf1(:,0), pdf2(:,0)).atx.x)&
               &             + (PartonLuminosity(grid, pdf1(:,0), P_x_pdf2(:,0)).atx.x))*cs%ln_muF2_M2
       end if
       if (cs%ln_muR2_M2 /= zero) then
          res(2) = res(2) + 2*twopi*beta0*cs%ln_muR2_M2 * res(1)
       end if


       res(2) = res(2) * cs%as2pi
       res    = res * cs%norm_fb
    case("DY")
       call InitDYCoeffs(grid)
       lumi_qqbar(0:grid%ny) = lumi_all_x(pdf1, pdf2, cs)
       res(1) = lumi_qqbar .atx. x

       res(2) = (DY_D1qg .conv. lumi_DY_qg(pdf1, pdf2, cs)) .atx. x
       res(2) = res(2) + ( (DY_D1qqbar * lumi_qqbar) .atx. x)

       if (cs%ln_muF2_M2 /= zero) then
          P_x_pdf1 = dglap_h%P_LO * pdf1
          P_x_pdf2 = dglap_h%P_LO * pdf2
          res(2) = res(2) - (  (lumi_all_x(P_x_pdf1, pdf2, cs).atx.x)&
               &             + (lumi_all_x(pdf1, P_x_pdf2, cs).atx.x))*cs%ln_muF2_M2
       end if
       ! NB: do nothing for muR scale dependence at this order.

       res(2) = res(2) * cs%as2pi

    end select

  end function cross_sections


  !----------------------------------------------------------------------
  function lumi_DY_qg(pdf1,pdf2,cs) result(res)
    real(dp),                     intent(in) :: pdf1(:,-6:), pdf2(:,-6:)
    type(process_and_parameters), intent(in) :: cs
    real(dp)                     :: res(size(pdf1,1)), weighted_q
    real(dp)                     :: weighted_q1(size(pdf1,1)), weighted_q2(size(pdf1,1))
    !----------------------------
    integer i
    if (cs%proc /= "DY") call wae_error("lumi_DY_qg: should not be called with proc != DY")

    weighted_q1 = zero
    weighted_q2 = zero
    do i = 1, nf_int
       weighted_q1 = weighted_q1 + gv2_ga2(i) * (pdf1(:,-i) + pdf1(:,i))
       weighted_q2 = weighted_q2 + gv2_ga2(i) * (pdf2(:,-i) + pdf2(:,i))
    end do

    res = (  PartonLuminosity(grid, weighted_q1(:), pdf2(:,0)) &
         & + PartonLuminosity(grid, weighted_q2(:), pdf1(:,0)) ) * cs%norm_fb
  end function lumi_DY_qg


end module pdfs_tools


!======================================================================
!======================================================================
!===== the following module contains routines which perform the  ======
!===== expansion of the luminosity factor up to order O(alphas^3) =====
!======================================================================
!======================================================================

module pdf_expansions
  use hoppet_v1; use rad_tools
  use pdfs_tools
  implicit none

  private :: Delete_PDFE
  interface Delete
    module procedure Delete_PDFE, Delete_LumiE
  end interface

  !- type to hold the expansion of a PDF;
  type PDFExpansion
    ! pij is coefficient of (as2pi(muR))^i * t^j, where t = ln muF1^2/muF0^2 = -2L
    real(dp), pointer :: p(:,:)
    real(dp), pointer :: p11(:,:)
    real(dp), pointer :: p21(:,:)
    real(dp), pointer :: p22(:,:)
    real(dp), pointer :: p31(:,:)
    real(dp), pointer :: p32(:,:)
    real(dp), pointer :: p33(:,:)

    ! cij is part of coefficient of (as2pi(muR))^i * t^j, where t =
    ! ln muF1^2/muF0^2 = -2L that includes a factor of the coefficient 
    ! function
    real(dp), pointer :: c10(:,:)
    real(dp), pointer :: c21(:,:)
    real(dp), pointer :: c31(:,:)
    real(dp), pointer :: c32(:,:)
  end type PDFExpansion

  !--------------------------
  type LumiExpansion
    real(dp), pointer :: l00(:) 

    real(dp), pointer :: l10(:) 
    real(dp), pointer :: l11(:) 

    real(dp), pointer :: l20(:) 
    real(dp), pointer :: dl21(:)
    real(dp), pointer :: l21(:) 
    real(dp), pointer :: l22(:) 

    real(dp), pointer :: l30(:) 
    real(dp), pointer :: l31(:) 
    real(dp), pointer :: l32(:) 
    real(dp), pointer :: l33(:) 
  end type LumiExpansion

  !--------------------------
  type LumiExpansion_at_x
    real(dp) :: l00 

    real(dp) :: l10 
    real(dp) :: l11 

    real(dp) :: l20 
    real(dp) :: dl21 
    real(dp) :: l21 
    real(dp) :: l22 

    real(dp) :: l30 
    real(dp) :: l31 
    real(dp) :: l32 
    real(dp) :: l33 
  end type LumiExpansion_at_x

  

contains

  !---------------------------------------------------
  subroutine AllocPDFExpansion(pdfe)
    type(PDFExpansion), intent(inout) :: pdfe
    call AllocPDF(grid,pdfe%p  )
    call AllocPDF(grid,pdfe%p11)
    call AllocPDF(grid,pdfe%p21)
    call AllocPDF(grid,pdfe%p22)
    call AllocPDF(grid,pdfe%p31)
    call AllocPDF(grid,pdfe%p32)
    call AllocPDF(grid,pdfe%p33)

    call AllocPDF(grid,pdfe%c10)
    call AllocPDF(grid,pdfe%c21)
    call AllocPDF(grid,pdfe%c31)
    call AllocPDF(grid,pdfe%c32)

  end subroutine AllocPDFExpansion
  subroutine Delete_PDFE(pdfe)
    type(PDFExpansion), intent(inout) :: pdfe
    deallocate(pdfe%p  )
    deallocate(pdfe%p11)
    deallocate(pdfe%p21)
    deallocate(pdfe%p22)
    deallocate(pdfe%p31)
    deallocate(pdfe%p32)
    deallocate(pdfe%p33)

    deallocate(pdfe%c10)
    deallocate(pdfe%c21)
    deallocate(pdfe%c31)
    deallocate(pdfe%c32)
  end subroutine Delete_PDFE

  !---------------------------------------------------
  subroutine AllocLumiExpansion(lumie)
    type(LumiExpansion), intent(inout) :: lumie
    call AllocGridQuant(grid,lumie%l00) 
    call AllocGridQuant(grid,lumie%l10) 
    call AllocGridQuant(grid,lumie%l11) 
    call AllocGridQuant(grid,lumie%l20) 
    call AllocGridQuant(grid,lumie%dl21) 
    call AllocGridQuant(grid,lumie%l21) 
    call AllocGridQuant(grid,lumie%l22) 
    call AllocGridQuant(grid,lumie%l30) 
    call AllocGridQuant(grid,lumie%l31) 
    call AllocGridQuant(grid,lumie%l32) 
    call AllocGridQuant(grid,lumie%l33) 
  end subroutine AllocLumiExpansion

  subroutine Delete_LumiE(lumie)
    type(LumiExpansion), intent(inout) :: lumie
    deallocate(lumie%l00) 
    deallocate(lumie%l10) 
    deallocate(lumie%l11) 
    deallocate(lumie%l20) 
    deallocate(lumie%dl21) 
    deallocate(lumie%l21) 
    deallocate(lumie%l22) 
    deallocate(lumie%l30) 
    deallocate(lumie%l31) 
    deallocate(lumie%l32) 
    deallocate(lumie%l33) 
  end subroutine Delete_LumiE
  
  

  !----------------------------------------------------------------------
  ! given an input pdf that has already been set in pdfe%p (assuming 
  ! the whole of pdfe has been allocated), 
  subroutine set_pdf_expansion(pdfe,u,cs)
    type(PDFExpansion), intent(inout) :: pdfe
    real(dp),           intent(in)    :: u
    !----------------------------------
    real(dp) :: P0_pdf(size(pdfe%p,1),size(pdfe%p,2))
    real(dp) :: P0_P0_pdf(size(pdfe%p,1),size(pdfe%p,2)), P1_pdf(size(pdfe%p,1),size(pdfe%p,2))
    real(dp) :: P2_pdf(size(pdfe%p,1),size(pdfe%p,2))
    real(dp) :: P1_P0_pdf(size(pdfe%p,1),size(pdfe%p,2)), P0_P1_pdf(size(pdfe%p,1),size(pdfe%p,2))
    real(dp) :: P0_P0_P0_pdf(size(pdfe%p,1),size(pdfe%p,2))

    type(process_and_parameters), intent(in)  :: cs
    type(split_mat)                           :: C1

    !--- set up the convolutions of splitting fn * pdf
    P0_pdf = dglap_h%P_LO * pdfe%p
    P0_P0_pdf = dglap_h%P_LO * P0_pdf
    P0_P0_P0_pdf = dglap_h%P_LO * P0_P0_pdf
    
    P1_pdf = dglap_h%P_NLO * pdfe%p
    P1_P0_pdf = dglap_h%P_NLO * P0_pdf
    P0_P1_pdf = dglap_h%P_LO * P1_pdf
    
    P2_pdf = dglap_h%P_NNLO * pdfe%p

    !--- now work out the various expansions
    pdfe%p11 = P0_pdf
    
    pdfe%p22 = half*P0_P0_pdf - pi*beta0 * P0_pdf
    pdfe%p21 = P1_pdf - twopi*u*beta0 * P0_pdf

    pdfe%p33 = (P0_P0_P0_pdf - (6*Pi*beta0)*P0_P0_pdf + (8*Pi**2*beta0**2)*P0_pdf)/6._dp
    pdfe%p32 = (P0_P1_pdf + P1_P0_pdf)/2._dp - (2*Pi*u*beta0)*P0_P0_pdf -  &
         &    (2*Pi*beta0)*P1_pdf + (4*Pi**2*u*beta0**2- 2*Pi**2*beta1)*P0_pdf
    pdfe%p31 = P2_pdf - (4*Pi*u*beta0)*P1_pdf + (4*Pi**2*u*(u*beta0**2-beta1))*P0_pdf
    
    ! these pieces should add as2pi(muF1/muF0 * muR) * C1 * pdf_evolved
    ! as2pi(muF1/muF0 * muR) = as2pi(muR) - b0*t*as2pi^2 
    call InitSplitMat(C1, C1_matrix)
    call AddWithCoeff(C1, dglap_h%P_LO, -cs%ln_muF2_M2+cs%ln_Q2_M2)

    pdfe%c10 = C1*pdfe%p
    pdfe%c21 = C1*pdfe%p11 - twopi*beta0*pdfe%c10
    pdfe%c31 = C1*pdfe%p21 ! - twopi**2*beta1*pdfe%c10 # this commented piece would appear if the two-loop running were used instead
    pdfe%c32 = C1*pdfe%p22 - twopi*beta0*(C1*pdfe%p11) + twopi**2*beta0**2*pdfe%c10
  end subroutine set_pdf_expansion

  !----------------------------------------------------------------------

  subroutine set_lumi_expansion(pdf1e, pdf2e, lumie, include_c1, include_h1, include_c1sq, cs)
    type(PDFExpansion),           intent(in)  :: pdf1e, pdf2e
    type(LumiExpansion),          intent(out) :: lumie
    logical,                      intent(in)  :: include_c1, include_c1sq, include_h1
    type(process_and_parameters), intent(in)  :: cs

    lumie%l00 = lumi_all_x(pdf1e%p, pdf2e%p, cs)
    lumie%l11 = lumi_all_x(pdf1e%p11, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%p11,cs)
    lumie%l21 = lumi_all_x(pdf1e%p21, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%p21,cs)
    lumie%l22 = lumi_all_x(pdf1e%p11, pdf2e%p11,cs) + lumi_all_x(pdf1e%p22, pdf2e%p,cs) +&
            & lumi_all_x(pdf1e%p, pdf2e%p22,cs)
    lumie%l31 = lumi_all_x(pdf1e%p31, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%p31,cs) 
    lumie%l32 = lumi_all_x(pdf1e%p32, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%p32,cs) + &
            & lumi_all_x(pdf1e%p21, pdf2e%p11,cs) + lumi_all_x(pdf1e%p11, pdf2e%p21,cs) 
    lumie%l33 = lumi_all_x(pdf1e%p33, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%p33,cs) + &
            & lumi_all_x(pdf1e%p22, pdf2e%p11,cs) + lumi_all_x(pdf1e%p11, pdf2e%p22,cs) 
     

    lumie%l10  = zero
    lumie%l20  = zero
    lumie%dl21 = zero
    lumie%l30  = zero
    
    ! we could imagine making things more efficient by first summing common parts, like
    ! c21 and p21, etc. before carrying out some of the calculations below
  
    if (include_c1) then
       if (include_h1) then
          lumie%l10  = lumie%l10  + H(1)*lumie%l00
          lumie%l31  = lumie%l31  + H(1)*lumie%l21
          lumie%dl21 = lumie%dl21 + H(1)*lumie%l11
          lumie%l21  = lumie%l21  + H(1)*lumie%l11
          lumie%l32  = lumie%l32  + H(1)*lumie%l22
       endif
      lumie%l10  = lumie%l10  + lumi_all_x(pdf1e%c10, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c10,cs)
      lumie%dl21 = lumie%dl21 + lumi_all_x(pdf1e%c21, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c21,cs)&
           &                  + lumi_all_x(pdf1e%c10, pdf2e%p11,cs) + lumi_all_x(pdf1e%p11, pdf2e%c10,cs)
      lumie%l21  = lumie%l21  + lumi_all_x(pdf1e%c21, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c21,cs)&
           &                  + lumi_all_x(pdf1e%c10, pdf2e%p11,cs) + lumi_all_x(pdf1e%p11, pdf2e%c10,cs)
      lumie%l31  = lumie%l31  + lumi_all_x(pdf1e%c31, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c31,cs)&
           &                  + lumi_all_x(pdf1e%c10, pdf2e%p21,cs) + lumi_all_x(pdf1e%p21, pdf2e%c10,cs)

      lumie%l32  = lumie%l32  + lumi_all_x(pdf1e%c32, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c32,cs)&
           &                  + lumi_all_x(pdf1e%c21, pdf2e%p11,cs) + lumi_all_x(pdf1e%p11, pdf2e%c21,cs)&
           &                  + lumi_all_x(pdf1e%c10, pdf2e%p22,cs) + lumi_all_x(pdf1e%p22, pdf2e%c10,cs)

      if (include_c1sq) then
        lumie%l20 = lumi_all_x(pdf1e%c10, pdf2e%c10,cs)
        lumie%l31 = lumie%l31 + lumi_all_x(pdf1e%c21, pdf2e%c10,cs) + lumi_all_x(pdf1e%c10, pdf2e%c21,cs)
     ! include terms containing H(1)   
        if (include_h1) then
           lumie%l20 = lumie%l20 + H(1)*(lumi_all_x(pdf1e%c10, pdf2e%p,cs) + lumi_all_x(pdf1e%p, pdf2e%c10,cs))
           lumie%l31 = lumie%l31 + H(1)*(lumi_all_x(pdf1e%c21, pdf2e%p,cs) &
                & + lumi_all_x(pdf1e%p, pdf2e%c21,cs) + lumi_all_x(pdf1e%c10, pdf2e%p11,cs) &
                & + lumi_all_x(pdf1e%p11, pdf2e%c10,cs))
           lumie%l30 = lumie%l30 + H(1)*lumi_all_x(pdf1e%c10, pdf2e%c10,cs)
        endif
      end if
    end if

  end subroutine set_lumi_expansion
  

!----------------------------------------------------------------------
!--- subroutine that evaluates the luminosity at a given x (threshold)
!----------------------------------------------------------------------

  subroutine set_lumi_expansion_at_x(lumi, lumie,cs)
    type(LumiExpansion_at_x),    intent(out) :: lumi
    type(LumiExpansion),       intent(inout) :: lumie
    type(process_and_parameters), intent(in) :: cs

    lumi%l00  = lumie%l00 .atx. (cs%M2_rts2 .with. grid)
    lumi%l10  = lumie%l10 .atx. (cs%M2_rts2 .with. grid)
    lumi%l11  = lumie%l11 .atx. (cs%M2_rts2 .with. grid)
    lumi%l20  = lumie%l20 .atx. (cs%M2_rts2 .with. grid)
    lumi%dl21 = lumie%dl21 .atx. (cs%M2_rts2 .with. grid)
    lumi%l21  = lumie%l21 .atx. (cs%M2_rts2 .with. grid)
    lumi%l22  = lumie%l22 .atx. (cs%M2_rts2 .with. grid)
    lumi%l30  = lumie%l30 .atx. (cs%M2_rts2 .with. grid)
    lumi%l31  = lumie%l31 .atx. (cs%M2_rts2 .with. grid)
    lumi%l32  = lumie%l32 .atx. (cs%M2_rts2 .with. grid)
    lumi%l33  = lumie%l33 .atx. (cs%M2_rts2 .with. grid)

  end subroutine set_lumi_expansion_at_x


end module pdf_expansions
