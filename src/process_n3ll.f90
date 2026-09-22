! SU(3), nf=5 Higgs HEFT and DY coefficients, in a=alpha_s/(2*pi).
! B3, A4 scheme conversion and H2: RadISH reference, not its prime branch.
! Exact cusp: arXiv:2002.04617 Eq. (6), including quartic Casimirs.
module process_n3ll
  use types
  use consts_dp
  use qcd
  use ew_parameters, only: mt
  use rad_tools, only: process_and_parameters,A,B,H,H1,CC
  use rapidity_n3ll
  implicit none
  private
  public :: coefficients_n3ll, cusp4_su3
contains
  function cusp4_su3(proc,nfl) result(v)
    character(len=*), intent(in) :: proc
    real(dp), intent(in) :: nfl
    real(dp) :: v,cr,df,da
    real(dp), parameter :: ca=3._dp,cf=4._dp/3,z5=1.03692775514336992633_dp
    select case(trim(proc))
    case('H'); cr=ca; df=15._dp/16; da=135._dp/8
    case('DY'); cr=cf; df=5._dp/36; da=5._dp/2
    case default; error stop 'N3LL cusp: unsupported process'
    end select
    v=nfl**3*cr*(64._dp/27*zeta3-32._dp/81) &
      +nfl**2*ca*cr*(-224._dp/15*zeta2**2+2240._dp/27*zeta3-608._dp/81*zeta2+923._dp/81) &
      +nfl**2*cf*cr*(64._dp/5*zeta2**2-640._dp/9*zeta3+2392._dp/81) &
      +nfl*ca**2*cr*(2096._dp/9*z5+448._dp/3*zeta3*zeta2-352._dp/15*zeta2**2 &
        -23104._dp/27*zeta3+20320._dp/81*zeta2-24137._dp/81) &
      +nfl*ca*cf*cr*(160*z5-128*zeta3*zeta2-352._dp/5*zeta2**2+3712._dp/9*zeta3 &
        +440._dp/3*zeta2-34066._dp/81) &
      +nfl*cf**2*cr*(-320*z5+592._dp/3*zeta3+572._dp/9) &
      +nfl*df*(-1280._dp/3*z5-256._dp/3*zeta3+256*zeta2) &
      +da*(-384*zeta3**2-7936._dp/35*zeta2**3+3520._dp/3*z5+128._dp/3*zeta3-128*zeta2) &
      +ca**3*cr*(-16*zeta3**2-20032._dp/105*zeta2**3-3608._dp/9*z5-352._dp/3*zeta3*zeta2 &
        +3608._dp/5*zeta2**2+20944._dp/27*zeta3-88400._dp/81*zeta2+84278._dp/81)
    ! Convert the coefficient of (alpha_s/(4*pi))^4.
    v=v/16
  end function

  subroutine coefficients_n3ll(cs,an,bn,betas,hard,delta2,delta3,algorithm)
    type(process_and_parameters), intent(in) :: cs
    character(len=*), intent(in) :: algorithm
    real(dp), intent(out) :: an(4),bn(3),betas(0:3),hard(2),delta2,delta3
    real(dp) :: hp2,aspow,lq,lr,lqr,f2
    real(dp), parameter :: zeta4=pi**4/90,zeta5=1.03692775514336992633_dp
    if(nf_int/=5 .or. abs(ca_def-3._dp)>1e-12_dp .or. abs(cf_def-4._dp/3)>1e-12_dp) &
      error stop 'N3LL process coefficients require SU(3), nf=5'
    if(cs%loop_mass/='none' .or. cs%observable/='ptj') &
      error stop 'N3LL requires ptj and loop-mass none'
    an(1:3)=A; bn(1:2)=B
    select case(trim(cs%proc))
    case('H')
      aspow=2
      hp2=5359._dp/54+137._dp/6*log(cs%M**2/mt**2)+1679._dp/24*pi**2+37._dp/8*pi**4-499._dp/6*zeta3
       bn(3) = (-270*cf_def*nf*(9*cf_def + 11*nf) + &
     &    ca_def*nf*(5*nf*(529 + 60*pi**2 - 144*zeta3) - &
     &       27*cf_def*(-1205 + 30*pi**2 + 8*pi**4 - 720*zeta4)) + &
     &    ca_def**2*nf*(-53185 - 5990*pi**2 + 738*pi**4 + 12360*zeta2 + 72720*zeta3 - 8100*zeta4) + &
     &    ca_def**3*(-2871*pi**4 + pi**2*(30545 - 5400*zeta3) + &
     &       5*(34219 - 91188*zeta3 + 24*zeta2*(-799 + 594*zeta3) - 12474*zeta4 + 38880*zeta5)))/9720._dp

    case('DY')
      aspow=0
      hp2=-57433._dp/972+281._dp/162*pi**2+22._dp/27*pi**4+1178._dp/27*zeta3
       bn(3) = ((ca_def*cf_def*nf*(-47.51165980795611_dp + (5188*Pi**2)/243._dp + (44*Pi**4)/45._dp - (3856*zeta3)/27._dp))/2._dp - &
            &    cf_def**2*nf*(63.370370370370374_dp - (8*Pi**4)/45._dp - (304*zeta3)/9._dp) - &
            &    ca_def*cf_def*nf*(85.90672153635117_dp - (412*Pi**2)/243._dp + (2*Pi**4)/27._dp - (904*zeta3)/27._dp) - cf_def*nf**2*(-2.5459533607681757_dp - (32*zeta3)/9._dp) + &
            &    (cf_def*nf**2*(26.52400548696845_dp - (80*Pi**2)/27._dp - (64*zeta3)/27._dp))/4._dp + &
            &    (cf_def**2*nf*(218.74074074074073_dp - (52*Pi**2)/9._dp - (56*Pi**4)/27._dp + (1024*zeta3)/9._dp))/2._dp + &
            &    ca_def**2*cf_def*(-95.5727023319616_dp - (7163*Pi**2)/243._dp - (83*Pi**4)/45._dp + (7052*zeta3)/9._dp - (88*Pi**2*zeta3)/9._dp - 272*zeta5) + &
            &    ca_def*cf_def**2*(-75.5_dp + (410*Pi**2)/9._dp + (494*Pi**4)/135._dp - (1688*zeta3)/3._dp - (16*Pi**2*zeta3)/3._dp - 240*zeta5) - &
            &    ca_def**2*cf_def*(-407.4471879286694_dp + (3196*Pi**2)/243._dp + (77*Pi**4)/135._dp + (12328*zeta3)/27._dp - (88*Pi**2*zeta3)/9._dp - 192*zeta5) + &
            &    cf_def**3*(-29 - 6*Pi**2 - (16*Pi**4)/5._dp - 136*zeta3 + (32*Pi**2*zeta3)/3._dp + 480*zeta5))/8._dp

    case default
      error stop 'N3LL supports only H and DY'
    end select
    an(4)=cusp4_su3(cs%proc,nf) + (CC*pi*(-288*beta1*pi*(28*nf+ca_def*(-202+189*zeta3)) &
      +beta0*(2*ca_def*nf*(-31313+3708*zeta2+12204*zeta3-2430*zeta4) &
      +nf*(32*nf*(58+81*zeta3)+27*cf_def*(-1711+912*zeta3+432*zeta4)) &
      +ca_def**2*(297029-332856*zeta3+72*zeta2*(-799+594*zeta3)-37422*zeta4+139968*zeta5))))/972._dp
    betas(0:2)=[beta0,beta1,beta2]
    betas(3)=(149753._dp/6+3564*zeta3-(1078361._dp/162+6508._dp/27*zeta3)*nf &
      +(50065._dp/162+6472._dp/81*zeta3)*nf**2+1093._dp/729*nf**3)/(4*pi)**4
    delta2=gamma_jet2(cs%jet_radius,CC,nf)-gamma_rsv2(CC,nf)
    delta3=delta_gamma3(cs%jet_radius,CC,nf,algorithm)
    f2=delta2/(32*A(1))
    lq=cs%ln_Q2_M2; lr=cs%ln_muR2_M2; lqr=cs%ln_Q2_muR2
    hard(1)=H(1)
    hard(2)=hp2+(A(1)*B(1)/2+A(1)*pi*beta0/3)*(-lq)**3+A(1)**2*lq**4/8 &
      +((B(1)**2-A(2))/2+pi*beta0*(B(1)+A(1)*lqr))*lq**2 &
      +(-B(2)+2*pi*beta0*B(1)*lqr)*(-lq)+B(1)*lq*H1-A(1)*lq**2*H1/2 &
      +4*aspow*((1+aspow)*pi**2*beta0**2*lr**2/2+pi**2*beta1*lr) &
      +H1*2*(1+aspow)*pi*beta0*lr+2*aspow*pi*beta0*lr*(-A(1)*lq/2+B(1))*lq &
      +8*A(1)*pi*beta0*zeta3/3-4*A(1)*lq*f2
  end subroutine
end module
