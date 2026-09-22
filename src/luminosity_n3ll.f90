! Unprimed N3LL jet-veto luminosity. No helicity-flip G terms for ptj.
! Hard coefficients must be supplied in the same alpha_s/(2*pi) scheme.
module luminosity_n3ll
  use types
  use consts_dp
  use hoppet_v1
  use rad_tools, only: process_and_parameters
  use pdfs_tools, only: grid,dglap_h,get_pdfs,lumi_at_x
  use coefficient_functions_n3ll, only: InitCoeffMatrix
  use coefficient_functions_deltaptj, only: InitCoeffMatrixDeltaPtj
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: init_luminosity_n3ll, n3ll_lumi, n3ll_prefactor
  public :: n3ll_lumi_series, pdf_series_n3ll
  public :: n3ll_lumi_from_pdfs
  type(split_mat), save :: c1,c2
  logical, save :: ready=.false.
  real(dp), save :: radius_saved
  integer, save :: nf_saved
contains
  subroutine init_luminosity_n3ll(radius,data_directory)
    real(dp), intent(in) :: radius
    character(len=*), intent(in) :: data_directory
    type(split_mat) :: delta,g1
    ! grid and dglap_h must already be initialized by pdfs_tools.
    if(ready) then
      call Delete(c1); call Delete(c2)
      ready=.false.
    end if
    call InitCoeffMatrix(grid,c1,c2,g1)
    call InitCoeffMatrixDeltaPtj(grid,delta,radius,data_directory)
    call AddWithCoeff(c2,delta,one)
    call Delete(delta); call Delete(g1)
    radius_saved=radius; nf_saved=nf_int; ready=.true.
  end subroutine

  ! l1 is the sum of C1 on either leg; l2 includes C2 on either leg
  ! and the C1*C1 product. The explicit hard/collinear prefactor is
  ! truncated at order a^2; only running produces higher-order terms.
  function n3ll_prefactor(l0,l1,l2,a,lambda,b0,b1,hard1,hard2) result(value)
    real(dp), intent(in) :: l0,l1,l2,a,lambda,b0,b1,hard1,hard2
    real(dp) :: value,ar1,ar2,w
    if(.not.all(ieee_is_finite([l0,l1,l2,a,lambda,b0,b1,hard1,hard2]))) &
      error stop 'N3LL luminosity: nonfinite input'
    if(a<0 .or. lambda<0 .or. lambda>=half .or. b0<=0) &
      error stop 'N3LL luminosity: invalid coupling or lambda'
    w=one-two*lambda
    ar1=a/w
    ar2=ar1*(one-twopi*a/w*b1/b0*log(w))
    value=(one+a*hard1+a**2*hard2)*l0 &
      +(ar2+a*ar1*hard1)*l1+ar1**2*l2
  end function

  function n3ll_lumi(muF,lambda,cs,hard1,hard2) result(value)
    real(dp), intent(in) :: muF,lambda,hard1,hard2
    type(process_and_parameters), intent(in) :: cs
    real(dp) :: value
    real(dp) :: p(0:grid%ny,-6:7),q(0:grid%ny,-6:7)
    if(.not.ready) error stop 'N3LL luminosity: call initializer first'
    if(nf_int/=nf_saved) error stop 'N3LL luminosity: active flavour count changed'
    if(cs%observable/='ptj') error stop 'N3LL luminosity supports only ptj'
    if(abs(cs%jet_radius-radius_saved)>1e-12_dp) error stop 'N3LL luminosity: radius changed'
    if(.not.ieee_is_finite(muF) .or. muF<=0) error stop 'N3LL luminosity: invalid PDF scale'
    call get_pdfs(muF,cs%collider,p,q)
    value=n3ll_lumi_from_pdfs(p,q,lambda,cs,hard1,hard2)
  end function

  function n3ll_lumi_from_pdfs(p,q,lambda,cs,hard1,hard2) result(value)
    real(dp), intent(in) :: p(0:grid%ny,-6:7),q(0:grid%ny,-6:7),lambda,hard1,hard2
    type(process_and_parameters), intent(in) :: cs
    real(dp) :: value,l0,l1,l2
    real(dp) :: c1p(0:grid%ny,-6:7),c1q(0:grid%ny,-6:7)
    real(dp) :: c2p(0:grid%ny,-6:7),c2q(0:grid%ny,-6:7)
    if(.not.ready) error stop 'N3LL luminosity: call initializer first'
    call coefficient_actions(p,cs,c1p,c2p)
    call coefficient_actions(q,cs,c1q,c2q)
    l0=lumi_at_x(p,q,cs)
    l1=lumi_at_x(c1p,q,cs)+lumi_at_x(p,c1q,cs)
    l2=lumi_at_x(c2p,q,cs)+lumi_at_x(p,c2q,cs)+lumi_at_x(c1p,c1q,cs)
    value=n3ll_prefactor(l0,l1,l2,cs%as2pi,lambda,beta0,beta1,hard1,hard2)
  end function

  subroutine coefficient_actions(p,cs,out1,out2)
    real(dp), intent(in) :: p(0:grid%ny,-6:7)
    type(process_and_parameters), intent(in) :: cs
    real(dp), intent(out) :: out1(0:grid%ny,-6:7),out2(0:grid%ny,-6:7)
    real(dp) :: p0(0:grid%ny,-6:7),c1p(0:grid%ny,-6:7),d,lqr
    d=cs%ln_Q2_M2-cs%ln_muF2_M2; lqr=cs%ln_Q2_muR2
    p0=dglap_h%P_LO .conv. p
    c1p=c1 .conv. p
    out1=c1p+d*p0
    out2=(c2 .conv. p)+(pi*beta0*d**2-twopi*beta0*d*lqr)*p0 &
      +d*(dglap_h%P_NLO .conv. p)+half*d**2*(dglap_h%P_LO .conv. p0) &
      +d*(c1 .conv. p0)-twopi*beta0*lqr*c1p
  end subroutine

  function n3ll_lumi_series(L,cs,hard1,hard2) result(v)
    real(dp), intent(in) :: L,hard1,hard2
    type(process_and_parameters), intent(in) :: cs
    real(dp) :: v(0:3),l0(0:3),l1(0:2),l2(0:1),r1,r2
    real(dp) :: p(0:grid%ny,-6:7,0:3),q(0:grid%ny,-6:7,0:3)
    real(dp) :: cp(0:grid%ny,-6:7,0:2),cq(0:grid%ny,-6:7,0:2)
    real(dp) :: dp2(0:grid%ny,-6:7,0:2),dq2(0:grid%ny,-6:7,0:2)
    integer :: i,j,k
    if(.not.ready) error stop 'N3LL expansion: initialize coefficient matrices first'
    if(nf_int/=nf_saved .or. abs(cs%jet_radius-radius_saved)>1e-12_dp) &
      error stop 'N3LL expansion: nf or radius changed'
    call get_pdfs(cs%muF,cs%collider,p(:,:,0),q(:,:,0))
    call pdf_series_n3ll(p,L,cs)
    call pdf_series_n3ll(q,L,cs)
    do i=0,2
      call coefficient_actions(p(:,:,i),cs,cp(:,:,i),dp2(:,:,i))
      call coefficient_actions(q(:,:,i),cs,cq(:,:,i),dq2(:,:,i))
    end do
    l0=0; l1=0; l2=0
    do k=0,3
      do i=0,k
        j=k-i
        l0(k)=l0(k)+lumi_at_x(p(:,:,i),q(:,:,j),cs)
        if(k<=2) l1(k)=l1(k)+lumi_at_x(cp(:,:,i),q(:,:,j),cs)+lumi_at_x(p(:,:,i),cq(:,:,j),cs)
        if(k<=1) l2(k)=l2(k)+lumi_at_x(dp2(:,:,i),q(:,:,j),cs) &
          +lumi_at_x(p(:,:,i),dq2(:,:,j),cs)+lumi_at_x(cp(:,:,i),cq(:,:,j),cs)
      end do
    end do
    r1=4*pi*beta0*L; r2=8*pi**2*beta1*L
    v(0)=l0(0)
    v(1)=l0(1)+l1(0)+hard1*l0(0)
    v(2)=l0(2)+l1(1)+r1*l1(0)+hard1*l0(1)+l2(0)+hard2*l0(0)+hard1*l1(0)
    v(3)=l0(3)+l1(2)+r1*l1(1)+(r1**2+r2)*l1(0)+hard1*l0(2) &
      +l2(1)+2*r1*l2(0)+hard2*l0(1)+hard1*(l1(1)+r1*l1(0))
  end function

  subroutine pdf_series_n3ll(p,L,cs)
    real(dp), intent(inout) :: p(0:grid%ny,-6:7,0:3)
    real(dp), intent(in) :: L
    type(process_and_parameters), intent(in) :: cs
    real(dp) :: p0(0:grid%ny,-6:7),p00(0:grid%ny,-6:7),p000(0:grid%ny,-6:7)
    real(dp) :: p1(0:grid%ny,-6:7),p01(0:grid%ny,-6:7),p10(0:grid%ny,-6:7),p2(0:grid%ny,-6:7)
    real(dp) :: t,u,b0,b1
    t=-2*L; u=cs%ln_muF2_M2-cs%ln_muR2_M2
    b0=twopi*beta0; b1=twopi**2*beta1
    p0=dglap_h%P_LO .conv. p(:,:,0)
    p00=dglap_h%P_LO .conv. p0
    p000=dglap_h%P_LO .conv. p00
    p1=dglap_h%P_NLO .conv. p(:,:,0)
    p01=dglap_h%P_LO .conv. p1
    p10=dglap_h%P_NLO .conv. p0
    p2=dglap_h%P_NNLO .conv. p(:,:,0)
    p(:,:,1)=t*p0
    p(:,:,2)=t*(p1-b0*u*p0)+t**2*(p00-b0*p0)/2
    p(:,:,3)=t*(p2-2*b0*u*p1+(b0**2*u**2-b1*u)*p0) &
      +t**2*((p01+p10)/2-b0*u*p00-b0*p1+(b0**2*u-b1/2)*p0) &
      +t**3/6*(p000-3*b0*p00+2*b0**2*p0)
  end subroutine
end module
