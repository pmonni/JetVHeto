! Three-loop RAD from the authors' Zenodo Wolfram tables.
! All gamma coefficients multiply (alpha_s/(4*pi))**n.
module rapidity_n3ll
  use types
  use consts_dp
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: gamma_rsv2, gamma_rsv3, gamma_jet2, delta_gamma3, gamma_jet3
  public :: jet_radius_exponent, jet_radius_series
  include 'rad_grid.inc'
contains
  function gamma_rsv2(cr,nf) result(v)
    real(dp), intent(in) :: cr,nf
    real(dp) :: v
    v=448*cr*nf*0.5_dp/27+3*cr*(-1616._dp/27+56*zeta3)
  end function
  function gamma_rsv3(cr,nf) result(v)
    real(dp), intent(in) :: cr,nf
    real(dp) :: v
    real(dp), parameter :: z5=1.03692775514336992633_dp
    v=3*cr*nf*0.5_dp*(250504._dp/729-1648*pi**2/243+8*pi**4/27-2144*zeta3/9) &
      +4._dp/3*cr*nf*0.5_dp*(6844._dp/27-32*pi**4/45-1216*zeta3/9) &
      +cr*nf**2*0.25_dp*(-14848._dp/729-256*zeta3/27) &
      +9*cr*(-594058._dp/729+6392*pi**2/243+154*pi**4/135+28528*zeta3/27 &
      -176*pi**2*zeta3/9-384*z5)
  end function
  function gamma_jet2(r,cr,nf) result(v)
    real(dp), intent(in) :: r,cr,nf
    real(dp) :: v,l2,l4
    l2=log(2._dp); l4=log(4._dp)
    if (.not.ieee_is_finite(r) .or. r<=0) error stop 'RAD: radius must be positive'
    v=-16._dp/3*cr**2*pi**2*r**2+4*cr**2*r**4 &
      +cr*nf*0.5_dp*(2960._dp/27-32*pi**2/9+3071*r**2/1350-168401*r**4/1587600 &
      +7001023*r**6/762048000-5664846191._dp*r**8/8851949568000._dp &
      -976*l2/9-56*r**2*l2/45+106*r**4*l2/945-11*r**6*l2/1575 &
      +4001*r**8*l2/7484400+232*l4/9-64*l2*l4/3+16*(-23+24*l2)*log(r)/9) &
      +3*cr*(-8056._dp/27+88*pi**2/9+1429*r**2/2700+4*pi**2*r**2/3 &
      -9383279*r**4/6350400-pi**2*r**4/54+74801417*r**6/1524096000 &
      -50937246539._dp*r**8/35407798272000._dp-pi**2*r**8/388800 &
      +2216*l2/9+208*r**2*l2/45+587*r**4*l2/1890-23*r**6*l2/1050 &
      +28529*r**8*l2/29937600-548*l4/9+176*l2*l4/3 &
      -8*(-131+12*pi**2+132*l2)*log(r)/9+88*zeta3)
  end function
  function delta_gamma3(r,cr,nf,algorithm,error) result(v)
    real(dp), intent(in) :: r,cr,nf
    character(len=*), intent(in) :: algorithm
    real(dp), optional, intent(out) :: error
    real(dp) :: v,colors(7),weights(4),x(4),xx,rr,analytic,var
    integer :: alg,j,k,lo,near
    select case(trim(algorithm))
    case('antikt'); alg=1
    case('kt'); alg=2
    case('ca'); alg=3
    case default; error stop 'RAD: algorithm must be antikt, kt or ca'
    end select
    ! Grid radii were stored at single precision; allow that input roundoff only.
    if (.not.ieee_is_finite(r)) error stop 'RAD: nonfinite radius'
    if (r<rad_radii(1)-1e-7_dp .or. r>rad_radii(rad_grid_size)+1e-7_dp) &
      error stop 'RAD: radius outside authors grid; extrapolation forbidden'
    rr=max(rad_radii(1),min(r,rad_radii(rad_grid_size)))
    near=minloc(abs(rad_radii-rr),dim=1)
    if (abs(rad_radii(near)-rr)<1e-7_dp) rr=rad_radii(near)
    lo=max(1,min(rad_grid_size-3,count(rad_radii<rr)-1))
    x=log(rad_radii(lo:lo+3)); xx=log(rr); weights=1
    do j=1,4
      do k=1,4
        if (k/=j) weights(j)=weights(j)*(xx-x(k))/(x(j)-x(k))
      end do
    end do
    colors=[9*cr,3*cr**2,cr**3,1.5_dp*cr*nf,2._dp/3*cr*nf, &
      0.5_dp*cr**2*nf,0.25_dp*cr*nf**2]
    v=0; var=0
    do j=1,4
      v=v+weights(j)*dot_product(rad_values(:,lo+j-1,alg),colors)
      ! Conservative fully correlated interpolation of the statistical errors.
      var=var+abs(weights(j))*sqrt(dot_product(rad_variances(:,lo+j-1,alg),colors**2))
    end do
    analytic=3*cr**2*(-2144*pi**2*r**2/27+32*pi**4*r**2/9+1064*r**4/9 &
      -8*pi**2*r**4/3-44*r**6/27-11*r**8/405-11*r**10/11340 &
      -11*r**12/243000-1408*r**2*zeta3/3) &
      +cr**2*nf*0.5_dp*(640*pi**2*r**2/27-352*r**4/9+16*r**6/27 &
      +4*r**8/405+r**10/2835+r**12/60750+512*r**2*zeta3/3)
    v=v+analytic
    if(present(error)) error=var
  end function
  function gamma_jet3(r,cr,nf,algorithm) result(v)
    real(dp), intent(in) :: r,cr,nf
    character(len=*), intent(in) :: algorithm
    real(dp) :: v
    v=gamma_rsv3(cr,nf)+delta_gamma3(r,cr,nf,algorithm)
  end function
  function jet_radius_exponent(L,a,betas,lqm,lqr,d2,d3,truncate_at_n3ll) result(v)
    real(dp), intent(in) :: L,a,betas(0:3),lqm,lqr,d2,d3
    logical, intent(in), optional :: truncate_at_n3ll
    real(dp) :: v,w,lambda,running,lr,scale_term,t
    lambda=twopi*a*betas(0)*L; w=1-2*lambda
    if(w<=0) error stop 'N3LL jet-radius exponent at or beyond Landau pole'
    ! RadISH inclusive (no-index) exponentiated prescription, retaining its
    ! beyond-N3LL products. d2=32*A1*f2 and d3=128*A1*f3; use our current
    ! anomalous dimensions, not RadISH's older numerical coefficients.
    ! The constant -d2*lqm*a^2/8 remains in H2. Do not add a separate
    ! beta0*lqr term to the running coupling: scale_term includes it.
    t=-twopi*a*betas(1)/betas(0)*log(w)/w
    running=1+t
    lr=lqm-lqr
    scale_term=((2-3*lambda+2*lambda**2)*lqm-lr)
    if(present(truncate_at_n3ll)) then
      if(truncate_at_n3ll) then
        ! Strict fixed-lambda N3LL: drop t^2, the Delta3 running products,
        ! and scale_term*t. Keep H2_R=-Delta2*lqm/8 in the hard function.
        v=d2/4*a**2*L/w**2*(1+2*t) &
          +d3/8*a**3*L/w**3 &
          -pi*betas(0)*d2*a**3*L/w**3*scale_term
        return
      endif
    endif
    v=d2/4*a**2*L/w**2*running**2 &
      +d3/8*a**3*L/w**3*running**3 &
      -pi*betas(0)*d2*a**3*L/w**3*scale_term*running
  end function

  function jet_radius_series(L,betas,lqm,lqr,d2,d3) result(v)
    real(dp), intent(in) :: L,betas(0:3),lqm,lqr,d2,d3
    real(dp) :: v(3)
    v(1)=0
    v(2)=d2*L/4
    v(3)=2*pi*betas(0)*d2*L**2-pi*betas(0)*d2*L*(lqm+lqr)+d3*L/8
  end function
end module rapidity_n3ll
