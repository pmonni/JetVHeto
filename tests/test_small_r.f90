program test_small_r
  use types
  use consts_dp
  use rad_tools
  use emsn_tools
  use pdfs_tools
  use resummation
  use expansion
  use matching
  use modified_logs
  use hoppet_v1, only: beta0
  implicit none
  type(process_and_parameters) :: cs,local,toy
  real(dp) :: pts(4),rs0(4),rs1(4),ex0(3,4),ex1(3,4),dl0(4),dl1(4),de0(2,4),de1(2,4)
  real(dp) :: ep(1),ll(1),lam(1),f(1),f2,f3,aa,err(2),born,original(1),al(1)
  real(dp) :: sig(0:3),bars(3,1),ex(3,1),dl(1),de(2,1),rs(1),matched(1)
  real(dp) :: s1,s2,s3,p1,p2,p3,c1,c2,c3,rem(2)
  integer :: ip,ir,i,k,n
  character(len=2048) :: pdf
  character(len=5) :: proc
  call get_command_argument(1,pdf)
  call get_command_argument(2,proc)
  if(len_trim(proc)==0) proc='H'
  call init_pdfs_from_LHAPDF(trim(pdf),0)
  call set_process_and_parameters(cs,'pp',trim(proc),13600._dp,125._dp,62.5_dp,62.5_dp,62.5_dp,5._dp,.4_dp,'none','ptj')
  cs%small_r=.false.;cs%small_r_ln2z=.false.;cs%small_r_R0=1
  cs%include_c1_squared=.false.;cs%matching_anew=.true.
  pts=[15._dp,30._dp,62.5_dp,125._dp]
  do ip=0,1
    cs%use_new_modlog=ip==1
    do ir=1,3
      cs%jet_radius=.2_dp*ir;cs%small_r_R0=.8_dp+.1_dp*ir
      call init_proc(cs)
      cs%small_r=.false.
      rs0=resummed_sigma(pts,cs,order_NNLL,dl0)
      ex0=expanded_sigma(pts,cs,order_NNLL,'anew',de0)
      cs%small_r=.true.
      rs1=resummed_sigma(pts,cs,order_NNLL,dl1)
      ex1=expanded_sigma(pts,cs,order_NNLL,'anew',de1)
      if(maxval(abs(dl1-dl0))>1e-12_dp) error stop 'small-R leaked into anew P'
      if(maxval(abs(de1-de0))>1e-12_dp) error stop 'small-R changed P1/P2'
      if(any(ex0(1:2,:)/=ex1(1:2,:))) error stop 'small-R changed expansion below a^3'
      do i=1,size(pts)
        if(cs%use_new_modlog) then
          call profile_point(cs,pts(i),local,ep)
        else
          local=cs;ep=pts(i);call init_proc(local)
        endif
        ll=Ltilde(ep/local%Q,local%p);lam=get_lambda(ll,local)
        f=small_r_factor(lam,local)
        al=local%as2pi/(1-2*lam)
        original=exp(-Rad_p(lam)*av_lnz_smallR(al,local%jet_radius,local%small_r_R0)) &
          +Rad_p(lam)*2*al*(non_incl(local%jet_radius,'all')-non_incl_lnR(local%jet_radius,local%small_r_R0))
        if(maxval(abs(f-original))>1e-14_dp) error stop 'Original small-R factor changed'
        original=Rad_p(lam)*2*al*non_incl(local%jet_radius,'all')
        if(abs(rs1(i)/rs0(i)-f(1)/exp(original(1)))>1e-12_dp) error stop 'anew radius factor routing'
        born=lumi_LL(local)
        f3=16*CC*non_incl_aslnR_sq(local%jet_radius,local%small_r_R0)*ll(1)
        if(abs(ex1(3,i)-ex0(3,i)-born*local%as2pi**3*f3)>1e-9_dp) error stop 'Cubic subtraction'
        if(ll(1)==0.and.f(1)/=1) error stop 'Small-R endpoint'
      enddo
    enddo
  enddo
  ! Independent numerical Taylor check of the actual all-order factor.
  call init_proc(cs)
  do ir=1,3
    toy=cs;toy%jet_radius=.2_dp*ir;toy%small_r_R0=1
    ll=1.3_dp
    f2=8*A(1)*non_incl(toy%jet_radius,'all')*ll(1)
    f3=8*pi*beta0*ll(1)*f2+16*CC*non_incl_aslnR_sq(toy%jet_radius,toy%small_r_R0)*ll(1)
    do k=1,2
      aa=.001_dp/2**(k-1);toy%as2pi=aa;toy%alphas_muR=twopi*aa
      lam=get_lambda(ll,toy);f=small_r_factor(lam,toy)
      err(k)=abs(f(1)-1-aa**2*f2-aa**3*f3)
    enddo
    if(err(1)/err(2)<14.or.err(1)/err(2)>18) error stop 'Factor remainder not O(a^4)'
    ! Matching algebra with the actual small-R factor and an independently
    ! specified Sudakov/luminosity, for each supported fixed-order accuracy.
    s1=-2;s2=-3;s3=-4;p1=.7_dp;p2=.4_dp;p3=.2_dp
    c1=s1+p1
    c2=s2+s1*s1/2+p2+s1*p1+f2
    c3=s3+s1*s2+s1**3/6+p3+s1*p2+(s2+s1*s1/2)*p1+f3+f2*c1
    do n=1,3
      do k=1,2
        aa=.001_dp/2**(k-1);toy%as2pi=aa;toy%alphas_muR=twopi*aa
        lam=get_lambda(ll,toy);f=small_r_factor(lam,toy)
        rs=f*exp(aa*s1+aa**2*s2+aa**3*s3)*(1+aa*p1+aa**2*p2+aa**3*p3)
        ex(:,1)=[aa*c1,aa**2*c2,aa**3*c3]
        dl=aa*p1+aa**2*p2+aa**3*p3;de(:,1)=[aa*p1,aa**2*p2]
        sig=[1._dp,aa*.3_dp,aa**2*.5_dp,aa**3*.8_dp];bars=0
        matched=match('anew',sig,ep,toy,rs,ex(:n,:),bars(:n,:),dl,de)
        rem(k)=abs(matched(1)-sum(sig(:n)))
      enddo
      if(rem(1)/rem(2)<.85_dp*2**(n+1).or.rem(1)/rem(2)>1.15_dp*2**(n+1)) &
        error stop 'Matched remainder has wrong perturbative order'
    enddo
  enddo
  toy%small_r_R0=toy%jet_radius
  f=small_r_factor(lam,toy)
  original=1+Rad_p(lam)*2*toy%as2pi/(1-2*lam)*non_incl(toy%jet_radius,'all')
  if(maxval(abs(f-original))>1e-14_dp) error stop 'R=R0 limit'
  print *, 'small-R: original factor, both logs, P separation, cubic expansion and NLO/NNLO/N3LO matching passed'
end program
