program jetvheto
  use sub_defs_io; use types; use consts_dp
  use rad_tools;  use pdfs_tools
  use expansion; use resummation; use reader
  use opts_cmdline_cardfile
  use ew_parameters
  use mass_corr
  use matching
  use banner
  use warnings_and_errors
  use interpolation
  implicit none

  ! --------------------------------------------------------------------
  type(process_and_parameters)  :: cs
  type(reader_opts)             :: opts, dummy_opts
  real(dp), allocatable, target :: pt_and_sigmabar(:,:), resummed(:), &
       &resummed_expanded(:,:), matched(:), interp_nnlo(:), pt_and_sigmabar_ref(:,:)
  real(dp), allocatable :: dlumi_lumi(:), dlumi_lumi_expanded(:,:)
  real(dp), allocatable :: final_result(:,:), splb(:), splc(:),spld(:)
  real(dp) :: jet_radius, ptmin, ptmax, ptcut, M, Q, muR, muF, p, rts, xQ, small_r_R0
  real(dp) :: sigma(0:3), norm_resum(1), sigma_native(0:1)
  real(dp), pointer  :: pt(:)
  integer            :: pdf_set, ibin, nbins, order, iunit = 7
  integer            :: fixed_order
  real(dp)           :: dy_error
  character(len=100) :: pdf_name
  character(len=2)   :: proc
  character(len=5)   :: collider
  character(len=20)  :: algo
  character(len=10)  :: matching_scheme
  character(len=200) :: outfile, infile, interp_infile
  character(len=4)   :: loop_mass
  character(len=3)   :: observable
  logical :: differential, check_logs, small_r, logscale, fivecol
  logical :: resum_only, do_matching, cross_section, onejet_cross_section
  logical :: give_error = .true.

  !---------------------------------------------------
  !! banners and help
  write(6,*)
  call print_banner
  write(6,*)
  call sleep(1)
  if (log_val_opt("-h") .or. log_val_opt("-help")) then
     call print_help
     stop
  end if

  !----------------------------------------------------------------------
  !! options dependent on whether doing matching, or just a resummation
  resum_only  = log_val_opt("-resum-only")
  ! add the following option to get the resummed cross section instead of the efficiency.
  ! valid only with matching scheme a and without the -resum-only option
  cross_section = log_val_opt("-cross-section",.false.)

  ! add the following option to get the inclusive 1 jet cross section
  onejet_cross_section = log_val_opt("-1j",.false.)

  do_matching = .not. resum_only
  if (do_matching) then
     infile = string_val_opt("-in", "")
     if (infile == '') then 
        write(0,*) 'Please specify fixed order file name with -in commandline option' 
        stop 
     endif
     ptmax = dble_val_opt("-ptmax",1e200_dp)
     call read_file_into_array(infile, pt_and_sigmabar, opts, max_col1 = ptmax)
     nbins = size(pt_and_sigmabar,2)

     ! total cross sections in nb
     sigma(0) = dble_val(opts,"xsct_lo")
     sigma(1) = dble_val(opts,"delta_xsct_nlo")
     sigma(2) = dble_val(opts,"delta_xsct_nnlo")
     sigma(3) = dble_val(opts,"delta_xsct_nnnlo",default=zero)

   ! internal use only: interpolate the NNLO column using some auxiliary file
     interp_infile = string_val_opt("-interp", "-1")
     if (interp_infile.ne."-1") then
        call read_file_into_array(interp_infile, pt_and_sigmabar_ref, dummy_opts, max_col1 = ptmax)
        ! total cross sections in nb
        allocate(interp_nnlo(nbins),splb(nbins),splc(nbins),spld(nbins))
        write(*,*) "interpolating NNLO from reference file ..."
        do ibin=1,nbins
           ! interpolate using cubic spline
!           call spline(pt_and_sigmabar_ref(1,:), pt_and_sigmabar_ref(3,:),splb,splc,spld,nbins)
!           interp_nnlo(ibin) = ispline(pt_and_sigmabar(1,ibin),pt_and_sigmabar_ref(1,:), &
!                & pt_and_sigmabar_ref(3,:),splb,splc,spld,nbins)
           ! interpolate using CERNLIB's polint
           call interpolate(pt_and_sigmabar_ref(1,:), pt_and_sigmabar_ref(3,:), &
                & pt_and_sigmabar(1,ibin), interp_nnlo(ibin), dy_error, 10)
        end do
     end if

  else
     ! set binning
     ptmin = dble_val_opt('-ptmin',1._dp)
     ptmax = dble_val_opt('-ptmax',300._dp)
     nbins = int_val_opt('-nbins', 300);

     ! create pt_and_sigmabar with zero entries for the sigmabar
     allocate(pt_and_sigmabar(3,nbins))
     pt => pt_and_sigmabar(1,:)
     if (log_val_opt('-log')) then
        do ibin = 1, nbins
           pt(ibin) = exp(log(ptmin)+(ibin-1)*(log(ptmax)-log(ptmin))/(nbins-1))
        end do
     else
        do ibin = 1, nbins
           pt(ibin) = ptmin+(ibin-1)*(ptmax-ptmin)/(nbins-1)
        end do
     end if
     pt_and_sigmabar(2:,:)=zero
     sigma = zero
  end if



  !----------------------------------------------------------------------
  ! set up the collider, process, PDF set and mass scale 
  collider = trim(string_val_opts_or_cmdline(opts,"collider","pp",give_error=give_error)) ! 'pp' or 'ppbar'
  rts      = dble_val_opts_or_cmdline(opts,"sqrts",13000._dp,give_error=give_error) ! Set the C.M. energy
  proc     = trim(string_val_opts_or_cmdline(opts,"proc","H",give_error=give_error)) ! 'H' or 'DY'
  if (proc == "H") then
     M = dble_val_opts_or_cmdline(opts,"M",125._dp,give_error=give_error) 
  else 
     M = MZ
  end if
  pdf_name = trim(string_val_opts_or_cmdline(opts,"pdf_name", &
       &                                     "NNPDF23_nnlo_as_0118.LHgrid",give_error=give_error))
  pdf_set  = int_val_opts_or_cmdline(opts,"pdf_set", 0, give_error=give_error)

  ! heavy-quark mass effects in Higgs production's virtual corrections (heavy top limit by default)
  loop_mass = trim(string_val_opts_or_cmdline(opts,"loop-mass","none",give_error=give_error)) ! 't', 't+b' or 'none' (large-mt limit)


  !--------- initialise PDFs -----------------------
  call init_pdfs_from_LHAPDF(pdf_name, pdf_set)


  !----------------------------------------------------------------------
  ! renormalisation, factorisation and resummation scales and other details of the resummation
  Q   = var_or_xVar("Q"  , xdefault=0.5_dp, opts=opts, M=M, override_allowed=.true.)
  muR = var_or_xVar("muR", xdefault=0.5_dp, opts=opts, M=M, override_allowed=.false.) ! must be consistent with fixed-order
  muF = var_or_xVar("muF", xdefault=0.5_dp, opts=opts, M=M, override_allowed=.false.) ! must be consistent with fixed-order

  ! the power for the damping of the logarithm at large pt (see Eq.24 of 1206.4998)
  p = dble_val_opts_or_cmdline(opts,"p",5._dp);
  ! the fixed order (1=NLO, 2=NNLO, 3=N3LO)
  fixed_order = int_val_opt('-fixed_order', fixed_order_NNLO)

  ! the resummation order (0 = LL, 1 = NLL, 2 = NNLL)
  order = int_val_opt("-order",order_NNLL)
  ! the matching scheme
  matching_scheme = string_val_opt("-scheme","a")

  ! the jet radius is that for any member of the longitudinally
  ! invariant generalised kt family (in so-called "inclusive" mode)
  ! with E-scheme recombination.
  jet_radius = dble_val_opts_or_cmdline(opts,'R',0.4_dp,give_error=give_error)
  ! the following line is for playing dangerous tests with R values that do 
  ! not match those in the fixed-order file.
  jet_radius = dble_val_opt("-badR",jet_radius)


  !----------------------------------------------------------------------
  ! decide where to send output and in what format
  outfile = string_val_opt("-out","/dev/stdout")
  if (outfile == "-") then
     iunit = 6
  else
     open(7,file=outfile)
  end if
  differential = log_val_opt('-differential',.false.) 
  ! for testing purposes only; NB: to effectively check the
  ! resummation's logarithmic structure requires long runs with 
  ! a program such as MCFM and in-depth understanding of the role
  ! of all the cutoffs that are present
  check_logs = log_val_opt("-check-logs",.false.)

  ! observable: default is pt_jet; other observable available is pt_boson (ptB)
  observable = string_val_opt("-observable","ptj") 
  if (observable == 'ptB' .and. fixed_order > 2) &
       &call wae_error('jetvheto: fixed order > 2 for ptB not supported') 

  ! the following cut is used as a lower cutoff for the resummed predictions
  ! its default value is MB/4 for ptB and 0 for ptj
  if (observable == 'ptB') then
     ptcut = dble_val_opt('-ptcut',M/four)
  else
     ptcut = dble_val_opt('-ptcut',zero)
  end if

  !! The following call is relevant for certain internal tests.
  !! Normal users should not need it.
  ! call small_alpha_init()

  !----------------------------------------------------------------------
  ! now initialise "cs" which is the object that contains the coefficients
  ! and is used also to communicate information on the collider, 
  ! scales, etc., throughout the rest of the program
  call set_process_and_parameters(cs, collider, proc, rts, M, muR, muF, Q, p, &
       &jet_radius,loop_mass,observable)
  ! this call arranges for the A and B resummation coefficients to be set
  call init_proc(cs) 
  ! final updates on parameters of cs
  cs%include_c1_squared = log_val_opt("-c1-squared",cs%include_c1_squared)

  ! - Added by FD -
  ! for small values of R, take into account higher orders of ln R with
  ! the -small-r option, setting small_r to true.
  ! small_r appears in resummed_sigma and expanded_sigma
  cs%small_r = log_val_opt("-small-r",.false.)
  ! Initial R0 for ln R resummation, set to 1 by default
  cs%small_r_R0 = dble_val_opt("-R0",1.0_dp)

  ! whether to include subleading ln^2z contribution
  cs%small_r_ln2z = log_val_opt("-small-r-ln2z",.false.)
  if (fixed_order >= fixed_order_N3LO .and. cs%small_r_ln2z) then
     call wae_error("fixed_order >= fixed_order_N3LO not yet compatible with -small-r-ln2z")
  end if
  
  
  !----------------------------------------------------------------------
  ! verify we've used all the command-line arguments and all the options
  ! from the fixed-order file
  if (.not. CheckAllArgsUsed(0)) stop
  call assert_all_opts_used(opts)

  !----------------------------------------------------------------------
  ! send documentation to the output files
  write(iunit,'(a)') "# "//trim(command_line())
  call print_parameters(iunit,cs)
  write(iunit,'(a)') "# order = "//order_string(order)
  write(iunit,'(a,i4)') "# pdf name and set = "//trim(pdf_name)//", ", pdf_set
  write(iunit,'(a)') "# matching scheme = "//trim(matching_scheme)

  !----------------------------------------------------------------------
  ! Get our own evaluation of the total cross section at LO and NLO.
  !
  ! The LO result will be used for ensuring consistent normalisations
  ! with external fixed order results; the NLO result could be used as a
  ! cross check of the user-supplied one, but this is not currently
  ! being done. Note that if finite top and bottom mass corrections are
  ! used, the LO cross section is computed exactly, while the NLO 
  ! cross section is computed in the large-mt limit! In this case one
  ! should always use the NLO value from the external fixed order result.
  sigma_native = cross_sections(cs)

  ! write cross sections here: note different units of measure for
  ! sigma native [fb] and what is provided by the input files [nb] 
  if (resum_only) then
     write(iunit,'(a,4f12.6)') '# total cross sections [nb] at LO, NLOdiff, NNLOdiff, N3LOdiff = ',sigma_native/1e6_dp
  else
     write(iunit,'(a,4f12.6)') '# total cross sections [nb] at LO, NLOdiff, NNLOdiff, N3LOdiff = ',sigma
  end if
  
  ! internal use only: overwrite the NNLO column with the interpolated one
  if ((loop_mass.ne."none").and.(interp_infile.ne."-1")) then
     ! rescale the large-mt limit according to HNNLO normalization
     pt_and_sigmabar(3,:) = interp_nnlo * norm_mass(cs%M, mt, mb, 't')
  end if


  if (.not. check_logs) then
     !-- perform the resummation, matching, and output of the results
     allocate(final_result(4,nbins))
     allocate(resummed(nbins), resummed_expanded(3,nbins), matched(nbins),&
          &dlumi_lumi(nbins),dlumi_lumi_expanded(2,nbins))

     ! the following calls get the resummation and its expansion, and
     ! also set dlumi_lumi and dlumi_lumi_expanded, which are
     ! ingredients for the matching below
     resummed          = resummed_sigma(pt_and_sigmabar(1,:), cs, order, dlumi_lumi)
     resummed_expanded = expanded_sigma(pt_and_sigmabar(1,:), cs, order, matching_scheme, dlumi_lumi_expanded)

     norm_resum =  sigma_native(0) 

     if (do_matching) then 
        resummed = resummed * sigma(0)/norm_resum(1)
        resummed_expanded = resummed_expanded * sigma(0)/norm_resum(1)
        matched = match(matching_scheme,&
             &  sigma, pt_and_sigmabar(1,:), cs, resummed, resummed_expanded(:fixed_order,:), &
             &  pt_and_sigmabar(2:fixed_order+1,:),dlumi_lumi,dlumi_lumi_expanded)

        ! -- Normalize efficiency for schemes a and b 
        if (matching_scheme == 'a') then 
           matched = matched / sum(sigma(0:fixed_order)) 
        elseif (matching_scheme == 'b') then 
           matched = matched / sum(sigma(0:fixed_order-1)) 
        elseif (matching_scheme == 'c') then
           matched = matched / sum(sigma(0:fixed_order-2))
        elseif (matching_scheme == 'cunx') then
           matched = matched / sum(sigma(0:fixed_order-2))
        elseif (matching_scheme == 'd') then
           matched = matched / sum(sigma(0:fixed_order-3))
        elseif (matching_scheme == 'moda') then 
           matched = matched / sum(sigma(0:fixed_order)) 
        elseif (matching_scheme == 'modRa') then 
           matched = matched / sum(sigma(0:fixed_order)) 
        elseif (matching_scheme == 'R') then
           matched = matched / sum(sigma(0:fixed_order)) 
        elseif (matching_scheme == 'loga') then
           matched = matched / matched(nbins)
        endif
     else
        sigma(0:1) = sigma_native(0:1)
        matched = zero 
     end if
     
     ! pt
     final_result(1,:) = pt_and_sigmabar(1,:)
     ! matched
     final_result(2,:) = matched
     ! pure resummed
     final_result(3,:) = resummed(:)/(sigma(0) + sum(sigma(1:order-1))) ! normalization  
     ! pure fixed order (in a scheme analogous to that used for the matching).
     final_result(4,:) = fixed_order_schemes(matching_scheme, sigma, pt_and_sigmabar(2:fixed_order+1,:))

     if (cross_section.and.do_matching) then
        if (matching_scheme == 'a') then
           final_result(2:,:) = final_result(2:,:) * sum(sigma(0:fixed_order)) 
           final_result(3,:)  = final_result(3,:)  * (sigma(0) + sum(sigma(1:order-1)))/sum(sigma(0:fixed_order))
           ! inclusive one-jet cross section
           if (onejet_cross_section) then
              final_result(2:,:) = sum(sigma(0:fixed_order)) - final_result(2:,:)
           end if
        else
           call wae_error('jetvheto: cross_section option is compatible only with matching_scheme a') 
        end if
     endif
        
     ! put out a header for the numbers
     write(iunit,'(a)') '#----------------------------------------------------------------------'
     if (differential) then
        if (cross_section.and.do_matching) then
           write(iunit,'(a)') '# the lines that follow give the derivative of the jet veto cross section wrt pt'
        else
           write(iunit,'(a)') '# the lines that follow give the derivative of the jet veto efficiencies wrt pt'
        end if
     else
        if (cross_section.and.do_matching) then
           write(iunit,'(a)') '# the lines that follow give the jet-veto cross section as a function of pt'
        else
           write(iunit,'(a)') '# the lines that follow give the jet-veto efficiencies as a function of pt'
        end if
     end if
     if (cross_section.and.do_matching) then
        write(iunit,'(a)') '#         pt            matched xs      resummed xs     fixed-order xs'
     else
        write(iunit,'(a)') '#         pt            matched eff      resummed eff     fixed-order eff'
     end if
  else
     !-- perform the resummation and output of results needed for the checks of
     !   logarithms
     allocate(final_result(7,nbins))
     allocate(resummed_expanded(3,nbins))

     ! col1 = pt, cols 2:3 are NLO and NNLO fixed-order integrated cross sections
     final_result(1:3,:) = pt_and_sigmabar(1:3,:)
     ! add in sigma(1) coefficient to NLO 
     final_result(2,:) = final_result(2,:) + sigma(1)

     ! get resummed expansion with right normalisation 
     resummed_expanded = expanded_sigma(pt_and_sigmabar(1,:), cs, order, matching_scheme)
     resummed_expanded = resummed_expanded * sigma(0)/sigma_native(0)
     ! subtract it from the fixed-order results
     final_result(2:3,:) = final_result(2:3,:) - resummed_expanded(1:2,:)
     ! and also store it in columns 4 & 5
     final_result(4,:) = resummed_expanded(1,:)
     final_result(5,:) = resummed_expanded(2,:)
     final_result(6,:) = resummed_expanded(3,:)
     final_result(7,:) = resummed_sigma(pt_and_sigmabar(1,:), cs, order, dlumi_lumi) 
     ! normalise all columns to sigma(0)
     final_result(2:,:) = final_result(2:,:) / sigma(0)
     write(iunit,'(a)') '# pt, LO(FO-resummed), NLO(FO-resummed), LO(exp. resum), NLO(exp. resum), NNLO(exp. resum), resummed'

  end if

  !---- final output of the results
  if (differential) then
     final_result(2:,:nbins-1) = (final_result(2:,2:nbins)-final_result(2:,1:nbins-1))/&
          &             spread(final_result(1,2:nbins)-final_result(1,1:nbins-1),ncopies=size(final_result,dim=1)-1,dim=1)
     final_result(1,:nbins-1) = (final_result(1,2:nbins)+final_result(1,1:nbins-1))/2
     do ibin=1,nbins-1
        if (final_result(1,ibin).lt.ptcut) cycle 
        !write(iunit,'(7es27.17e3)') final_result(:,ibin)
        write(iunit,'(5es18.8e3)') final_result(:,ibin)
     end do
  else
     do ibin=1,nbins
        if (final_result(1,ibin).lt.ptcut) cycle 
        !write(iunit,'(7es27.17e3)') final_result(:,ibin)
        write(iunit,'(5es18.8e3)') final_result(:,ibin)
     end do
  end if

  if (iunit /= 6) close(iunit)

contains
  

  !----------------------------------------------------------------------
  !! The following subroutine is relevant for certain internal tests.
  !! Normal users should not need it.
  subroutine small_alpha_init()
    ! Variables to check small alphas expansion 
    real(dp) :: alphas_new, alphas_ratio, alphasPDF
    
    ! Useful quantities for checks of the small alphas expansion
    alphas_new = dble_val_opt("-alphas_new",alphasPDF(muF))
    alphas_ratio = alphas_new/alphasPDF(muF)
    call init_pdfs_hoppet_evolution(pdf_name, pdf_set, muF, alphas_new,rts)
    !call assert_eq(muR, muF)
    sigma(1:2) = sigma(1:2) * (/ alphas_ratio, alphas_ratio**2 /)
    pt_and_sigmabar(2,:) = pt_and_sigmabar(2,:)*alphas_ratio
    pt_and_sigmabar(3,:) = pt_and_sigmabar(3,:)*alphas_ratio**2
    write(*,*) " new alphas(muF) = ", alphasPDF(muF)
  end subroutine small_alpha_init

end program jetvheto
