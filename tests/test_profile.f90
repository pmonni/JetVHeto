program test_profile
  use types
  use consts_dp
  use rad_tools
  use pdfs_tools
  use jetveto_n3ll
  use matching
  implicit none
  type(process_and_parameters) :: cs,local
  real(dp) :: pt(6),q,muM,xm,mass,x,ref,got,rs(6),ex(3,6),dl(6),de(2,6)
  real(dp) :: onept(1),r1(1),e1(3,1),d1(1),d2(2,1),hb(2),oldrs(6),oldex(3,6)
  real(dp) :: sig(0:3),bars(3,6),matched(6),eps
  integer :: i,j,k,n,mode
  character(len=2048) :: data,pdf
  call get_command_argument(1,data)
  call get_command_argument(2,pdf)
  mass=125
  do j=1,3
    xm=.25_dp*(j+1); muM=xm*mass
    do k=1,3
      q=mass*2._dp**(k-3)
      do i=0,200
        x=mass*i/100
        got=profile_scale(x,mass,q,xm)
        ref=reference(x,mass,q,muM)
        if(abs(got-ref)>1e-10_dp*mass) error stop 'Profile differs from RadISH polynomial'
        if(got<min(x,mass)-1e-12_dp) error stop 'Negative profiled log'
      enddo
      eps=1e-5_dp*mass
      if(abs(log(profile_scale(muM,mass,q,xm)/muM))>1e-12_dp) error stop 'Matching join'
      if(abs(log(profile_scale(mass,mass,q,xm)/mass))>1e-12_dp) error stop 'Endpoint'
      if(abs(log(profile_scale(mass-eps,mass,q,xm)/(mass-eps)))>1e-8_dp) error stop 'Endpoint slope'
    enddo
  enddo
  call init_pdfs_from_LHAPDF(trim(pdf),0)
  call set_process_and_parameters(cs,'pp','H',13600._dp,mass,62.5_dp,62.5_dp,62.5_dp,5._dp,.4_dp,'none','ptj')
  cs%small_r=.false.; cs%small_r_ln2z=.false.; cs%include_c1_squared=.false.
  call init_proc(cs)
  call initialize_jetveto_n3ll(cs,trim(data))
  pt=[5._dp,30._dp,62.5_dp,90._dp,125._dp,150._dp]
  do j=1,3
    cs%xM=.25_dp*(j+1); cs%use_new_modlog=.true.
    ! Independent cutoff checks below, at and above both profile joins.
    muM=cs%xM*mass
    do k=1,3
      cs%Q=mass*2._dp**(k-3)
      do i=1,7
        select case(i)
        case(1); x=muM-eps
        case(2); x=muM
        case(3); x=muM+eps
        case(4); x=(muM+mass)/2
        case(5); x=mass-eps
        case(6); x=mass
        case(7); x=mass+eps
        end select
        call profile_point(cs,x,local,onept)
        q=profile_scale(x,mass,cs%Q,cs%xM)
        if(local%Q/=q) error stop 'Cutoff changed coefficient scale'
        if(abs(local%ln_Q2_M2-2*log(q/mass))>1e-14_dp) error stop 'Cutoff changed scale log'
        got=log(local%Q/onept(1))
        if(x>muM) then
          if(got/=0._dp) error stop 'Nonzero logarithm above matching cutoff'
        else
          if(onept(1)/=x) error stop 'Cutoff changed sub-boundary evaluation point'
          if(abs(got-log(q/x))>1e-14_dp) error stop 'Incorrect logarithm below cutoff'
        endif
      enddo
    enddo
    cs%Q=62.5_dp
    call init_proc(cs)
    do mode=0,0
      cs%matching_anew=.true.
      hb=H
      rs=resum_jetveto_n3ll(pt,cs,dl)
      ex=expand_jetveto_n3ll(pt,cs,de)
      if(maxval(abs(H-hb))>1e-12_dp) error stop 'Profile leaked global coefficients'
      do i=1,size(pt)
        local=cs; local%use_new_modlog=.false.; local%p=-1
        local%Q=reference(pt(i),mass,cs%Q,cs%xM*mass)
        ! The unfactorized reference can round slightly below pt at its zero.
        local%Q=max(local%Q,min(pt(i),mass))
        local%ln_Q2_M2=2*log(local%Q/mass)
        local%ln_Q2_muR2=2*log(local%Q/local%muR)
        call init_proc(local)
        onept=min(pt(i),mass)
        if(pt(i)>cs%xM*mass) onept=local%Q
        r1=resum_jetveto_n3ll(onept,local,d1)
        e1=expand_jetveto_n3ll(onept,local,d2)
        if(abs(r1(1)-rs(i))>1e-7_dp.or.maxval(abs(e1(:,1)-ex(:,i)))>1e-7_dp) error stop 'Profile routing'
        if(abs(d1(1)-dl(i))>1e-10_dp.or.maxval(abs(d2(:,1)-de(:,i)))>1e-10_dp) error stop 'Profile prefactor'
      enddo
      call init_proc(cs)
      if(mode==0) then
        oldrs=rs; oldex=ex
      else
        if(maxval(abs(rs-oldrs))>1e-8_dp.or.maxval(abs(ex-oldex))>1e-8_dp) error stop 'Matching changed resummation'
      endif
      sig=[lumi_LL(cs),10000._dp,3000._dp,1000._dp]; bars=-20
      do n=1,3
        matched=match('a',sig,pt,cs,rs,ex(:n,:),bars(:n,:),dl,de)
        ! N3LL matched only to NLO retains its order-a^2 hard constant.
        if(n==1) matched(5:)=matched(5:)-ex(2,5:)
        if(maxval(abs(matched(5:)-sum(sig(:n))-sum(bars(:n,5))))>1e-7_dp) &
          error stop 'Profile matching does not recover fixed order at endpoint'
      enddo
    enddo
  enddo
  print *, '1809 profile polynomial checks; cutoff boundaries, routing, anew prefactor, endpoint checks passed'
contains
  function reference(v,a,q,b) result(y)
    real(dp),intent(in)::v,a,q,b
    real(dp)::y,m0,m1,m2,m3,den
    if(v>=a) then
      y=a
    else if(v<b) then
      y=q+3*(-4*q+2*b)/(2*b**2)*v**2+(8*q-3*b)/b**3*v**3-(3*q-b)/b**4*v**4
    else
      ! Original RadISH expressions with p=2, intentionally not factorized.
      den=-2*a*b*(a*a-b*b)+2*(a-b)*(a**3+b**3)
      m0=a*b*(a**3-3*a*a*b+3*a*b*b-b**3)/den
      m1=a*b*(-9*a*a+16*a*b-9*b*b+a**4/b**2+b**4/a**2)/((a-b)*den)
      m2=(-2*a**6*b+6*a**5*b**2-4*a**4*b**3-4*a**3*b**4+6*a**2*b**5-2*a*b**6) &
        /(a*a*(a-b)*b*b*den)
      m3=(a**5*b-4*a**4*b**2+6*a**3*b**3-4*a**2*b**4+a*b**5)/(a*a*(a-b)*b*b*den)
      y=m0+m1*v*v+m2*v**3+m3*v**4
    endif
  end function
end program
