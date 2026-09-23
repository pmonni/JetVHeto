program test_n3ll
  use types
  use consts_dp
  use hoppet_v1
  use rad_tools
  use pdfs_tools
  use process_n3ll
  use jetveto_n3ll
  use luminosity_n3ll
  use radiator_n3ll
  use matching
  implicit none
  type(process_and_parameters) :: cs
  real(dp), allocatable :: p(:,:),q(:,:),pstart(:,:),qstart(:,:)
  real(dp) :: an(4),bn(3),betas(0:3),hard(2),d2,d3,L,a0,lam,rad0(1),full,born
  real(dp) :: expanded(3,1),dl(2,1),errors(2),pt(1),pref,f,radiusc(3),rerr(2)
  real(dp) :: sig(0:3),bars(3,1),resum(1),dre(2,1),dres(1),matched(1),e1,e2,e3,p1,p2
  real(dp) :: scan(2,3),scan3(2,3),scale,w,runfac,amu,rp,sd,reference
  real(dp) :: expansion_reference(3,1),prefactor_reference(2,1)
  character(len=2048) :: directory,pdfname
  integer :: i,j,iproc,n
  call get_command_argument(1,directory)
  call get_command_argument(2,pdfname)
  if(len_trim(pdfname)==0) error stop 'Test requires a PDF set'
  call init_pdfs_from_LHAPDF(trim(pdfname),0)
  allocate(p(0:grid%ny,-6:7),q(0:grid%ny,-6:7),pstart(0:grid%ny,-6:7),qstart(0:grid%ny,-6:7))
  do iproc=1,2
    do j=-1,1
      call set_process_and_parameters(cs,'pp','H',13000._dp,125._dp,125._dp,125._dp, &
        125._dp,5._dp,0.4_dp,'none','ptj')
      if(iproc==2) then
        cs%proc='DY'; cs%norm_fb=1 ! Overall normalization irrelevant to relative-order test.
      end if
      ! Independently vary resummation, renormalization and factorization logs.
      cs%ln_Q2_M2=j*log(4._dp); cs%ln_Q2_muR2=-j*log(4._dp)
      cs%ln_muR2_M2=cs%ln_Q2_M2-cs%ln_Q2_muR2
      cs%ln_muF2_M2=0
      cs%small_r=.false.; cs%small_r_ln2z=.false.; cs%include_c1_squared=.false.
      call init_proc(cs)
      if(j==-1 .and. iproc==1) call initialize_jetveto_n3ll(cs,trim(directory))
      call coefficients_n3ll(cs,an,bn,betas,hard,d2,d3,'antikt')
      ! Independent RadISH no-index product form, at non-negligible coupling
      ! and with both signs of independent scale logs, for H and DY.
      do n=0,3
        L=real(n,dp); a0=0.02_dp
        lam=twopi*a0*betas(0)*L; w=1-2*lam
        runfac=1-twopi*a0*betas(1)/betas(0)*log(w)/w
        amu=a0/w*runfac; rp=4*an(1)*a0*L/w
        sd=-4*twopi*betas(0)*a0**2*d2/(32*an(1)) &
          *((-2+(3-2*lam)*lam)*(-cs%ln_Q2_M2)-cs%ln_muR2_M2)/w**2
        reference=rp*runfac*(2*amu*d2/(32*an(1))+4*amu**2*d3/(128*an(1))+sd)
        f=jet_radius_exponent(L,a0,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3)
        if(abs(f-reference)>1e-13_dp*max(1._dp,abs(reference))) error stop 'RadISH radius prescription mismatch'
        if(n==0.and.f/=0) error stop 'Radius endpoint not zero'
        ! Independent original strict expression, including scale compensation.
        reference=d2*a0**2*L/(4*w**2)*(1-4*pi*a0/w* &
          (betas(1)/betas(0)*log(w)+betas(0)*cs%ln_Q2_muR2)) &
          +d3*a0**3*L/(8*w**3)-d2*cs%ln_Q2_M2*a0**2/8*(1/w**2-1)
        f=jet_radius_exponent(L,a0,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3,.true.)
        if(abs(f-reference)>1e-13_dp*max(1._dp,abs(reference))) error stop 'Strict radius prescription mismatch'
        if(n==0.and.f/=0) error stop 'Strict radius endpoint not zero'
      enddo
      if(j==0) write(*,*) 'Process/cusp4/A4/B3/H2:',trim(cs%proc),cusp4_su3(cs%proc,5._dp),an(4),bn(3),hard(2)
      pt=[30._dp]; L=Ltilde_scalar(pt(1)/cs%Q,cs%p)
      call get_pdfs(cs%muF,cs%collider,pstart,qstart)
      born=lumi_at_x(pstart,qstart,cs)
      do i=1,2
        a0=0.0005_dp/2**(i-1); cs%as2pi=a0; cs%alphas_muR=twopi*a0
        expanded=expand_jetveto_n3ll(pt,cs,dl)
        ! Independent RK4 integration of DGLAP, not the polynomial PDF routine.
        p=pstart; q=qstart
        call evolve(p,a0,L,cs%ln_muF2_M2-cs%ln_muR2_M2)
        call evolve(q,a0,L,cs%ln_muF2_M2-cs%ln_muR2_M2)
        lam=twopi*a0*beta0*L
        pref=n3ll_lumi_from_pdfs(p,q,lam,cs,hard(1),hard(2))
        rad0=Rad([L],cs,order_NNLL)
        f=jet_radius_exponent(L,a0,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3)
        full=pref*exp(rad0(1)+4*a0**2*g4_n3ll(lam,an,bn,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2)+f)
        errors(i)=abs((full-born-sum(expanded))/born)
        radiusc=jet_radius_series(L,betas,cs%ln_Q2_M2,cs%ln_Q2_muR2,d2,d3)
        rerr(i)=abs(f-sum(radiusc*[a0,a0**2,a0**3]))
      end do
      write(*,*) 'Full expansion ratios ',trim(cs%proc),j,errors(1)/errors(2),rerr(1)/rerr(2)
      write(*,*) 'Relative expansion residuals:',errors
      if(abs(errors(1)/errors(2)-16._dp)>2) error stop 'Full alpha_s^3 expansion mismatch'
      if(abs(rerr(1)/rerr(2)-16._dp)>2) error stop 'Jet radius expansion mismatch'
    end do
  end do
  ! Both log choices: the runtime flag must not change any matching coefficient.
  do n=0,1
    cs%use_new_modlog=n==1
    cs%xM=0.5_dp
    do i=1,3
      pt=30._dp*i
      cs%truncate_rapidity_rge_at_n3ll=.false.
      expansion_reference=expand_jetveto_n3ll(pt,cs,prefactor_reference)
      cs%truncate_rapidity_rge_at_n3ll=.true.
      expanded=expand_jetveto_n3ll(pt,cs,dl)
      if(any(expanded/=expansion_reference).or.any(dl/=prefactor_reference)) &
        error stop 'Radius truncation changed matching expansion'
    enddo
  enddo
  cs%truncate_rapidity_rge_at_n3ll=.false.
  cs%use_new_modlog=.false.
  print *, 'Strict radius: independent formula and unchanged matching expansion checks passed'
  ! Q must cancel through order a^2 (up to negligible profile power terms).
  do iproc=1,2
    do j=1,3
      scale=125._dp*2._dp**(j-2)
      call set_process_and_parameters(cs,'pp','H',13000._dp,125._dp,125._dp,125._dp, &
        scale,20._dp,0.4_dp,'none','ptj')
      if(iproc==2) then
        cs%proc='DY'; cs%norm_fb=1
      end if
      call init_proc(cs)
      cs%as2pi=0.001_dp; cs%alphas_muR=twopi*cs%as2pi
      expanded=expand_jetveto_n3ll([5._dp],cs)
      scan(:,j)=expanded(1:2,1)/lumi_LL(cs)
      scan3(1,j)=expanded(3,1)/lumi_LL(cs)
      expanded=expand_jetveto_n3ll([10._dp],cs)
      scan3(2,j)=expanded(3,1)/lumi_LL(cs)
    end do
    write(*,*) 'Q cancellation through a^2 ',trim(cs%proc),maxval(abs(scan(:,1)-scan(:,2))), &
      maxval(abs(scan(:,3)-scan(:,2)))
    if(maxval(abs(scan(:,1)-scan(:,2)))>1e-9_dp .or. maxval(abs(scan(:,3)-scan(:,2)))>1e-9_dp) &
      error stop 'Resummation scale does not cancel through a^2'
    write(*,*) 'Q cancellation of cubic logarithms ',trim(cs%proc), &
      maxval(abs(scan3(1,:)-scan3(2,:)-(scan3(1,2)-scan3(2,2))))
    if(maxval(abs(scan3(1,:)-scan3(2,:)-(scan3(1,2)-scan3(2,2))))>1e-9_dp) &
      error stop 'Resummation scale does not cancel in cubic logarithms'
  end do
  ! Independent synthetic test of native matching scheme a through N3LO.
  ! Resum=e^(a e1+a^2 e2+a^3 e3)*(1+a p1+a^2 p2).
  e1=-3; e2=2; e3=7; p1=4; p2=9
  do n=1,3
    do i=1,2
      a0=0.002_dp/2**(i-1)
      sig=[1._dp,2*a0,5*a0**2,11*a0**3]; bars=0
      resum=[exp(a0*e1+a0*a0*e2+a0**3*e3)*(1+a0*p1+a0*a0*p2)]
      dres=[a0*p1+a0*a0*p2]; dre(:,1)=[a0*p1,a0*a0*p2]
      expanded(:,1)=[a0*(e1+p1),a0*a0*(e2+e1*e1/2+e1*p1+p2), &
        a0**3*(e3+e1*e2+e1**3/6+(e2+e1*e1/2)*p1+e1*p2)]
      matched=match('a',sig,pt,cs,resum,expanded(:n,:),bars(:n,:),dres,dre)
      errors(i)=abs(matched(1)-sum(sig(:n)))
    end do
    write(*,*) 'Matching order/remainder ratio:',n,errors(1)/errors(2)
    if(abs(errors(1)/errors(2)-2._dp**(n+1))>1) error stop 'Matching does not recover fixed order'
  end do
contains
  function Ltilde_scalar(x,power) result(y)
    real(dp), intent(in) :: x,power
    real(dp) :: y,tmp(1)
    tmp=Ltilde([x],power); y=tmp(1)
  end function
  subroutine evolve(pdf,a0,L,u)
    real(dp), intent(inout) :: pdf(0:grid%ny,-6:7)
    real(dp), intent(in) :: a0,L,u
    real(dp) :: k1(0:grid%ny,-6:7),k2(0:grid%ny,-6:7),k3(0:grid%ny,-6:7),k4(0:grid%ny,-6:7),h,t
    integer :: n
    h=-2*L/32
    do n=0,31
      t=n*h
      k1=rhs(pdf,a0,u+t)
      k2=rhs(pdf+h*k1/2,a0,u+t+h/2)
      k3=rhs(pdf+h*k2/2,a0,u+t+h/2)
      k4=rhs(pdf+h*k3,a0,u+t+h)
      pdf=pdf+h*(k1+2*k2+2*k3+k4)/6
    end do
  end subroutine
  function rhs(pdf,a0,t) result(v)
    real(dp), intent(in) :: pdf(0:grid%ny,-6:7),a0,t
    real(dp) :: v(0:grid%ny,-6:7),ar
    ar=a0-twopi*beta0*t*a0**2+(twopi**2*beta0**2*t**2-twopi**2*beta1*t)*a0**3
    v=ar*(dglap_h%P_LO .conv. pdf)+ar**2*(dglap_h%P_NLO .conv. pdf)+ar**3*(dglap_h%P_NNLO .conv. pdf)
  end function
end program
