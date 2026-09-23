program test_rad
  use types
  use rapidity_n3ll
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(dp) :: r,v,e,cr,nf,a1,legacy
  integer :: i,j,k,m
  character(len=20) :: mode
  character(len=6), parameter :: algs(3)=[character(len=6)::'antikt','kt','ca']
  real(dp), parameter :: radii(3)=[0.1_dp,0.4_dp,0.8_dp]
  real(dp), parameter :: expected(3)=[77165._dp,28537._dp,6645._dp]
  real(dp), parameter :: expected2(3)=[2174.69_dp,1135.60_dp,498.59_dp]
  call get_command_argument(1,mode)
  select case(trim(mode))
  case('low'); v=delta_gamma3(0.04_dp,3._dp,5._dp,'antikt')
  case('high'); v=delta_gamma3(1.01_dp,3._dp,5._dp,'antikt')
  case('nan'); v=delta_gamma3(ieee_value(0._dp,ieee_quiet_nan),3._dp,5._dp,'antikt')
  case('bad-algorithm'); v=delta_gamma3(0.4_dp,3._dp,5._dp,'unknown')
  case('dump')
    do j=1,3
      do k=1,2
        cr=3; if(k==2) cr=4._dp/3
        do m=1,3
          nf=0; if(m==2) nf=5; if(m==3) nf=6
          do i=1,20
            r=0.05_dp*i
            v=delta_gamma3(r,cr,nf,trim(algs(j)),e)
            write(*,'(a,1x,8(es25.16e3,1x))') trim(algs(j)),r,cr,nf,v,e, &
              gamma_jet2(r,cr,nf),gamma_rsv2(cr,nf),gamma_rsv3(cr,nf)
          end do
        end do
      end do
    end do
    stop
  case('')
  case default
    error stop 'Unknown test mode'
  end select
  if(len_trim(mode)>0) error stop 'Invalid input unexpectedly accepted'
  do i=1,3
    r=radii(i); v=gamma_jet3(r,3._dp,5._dp,'antikt')
    write(*,*) r,gamma_jet2(r,3._dp,5._dp),v
    if(abs(v-expected(i))>2._dp) error stop 'Published RAD benchmark mismatch'
    if(abs(gamma_jet2(r,3._dp,5._dp)-expected2(i))>0.005_dp) &
      error stop 'Published two-loop RAD benchmark mismatch'
  end do
  v=delta_gamma3(0.425_dp,3._dp,5._dp,'kt',e)
  if(e<=0 .or. abs(v)>1e6_dp) error stop 'Invalid interpolated result'
  if(abs(v-delta_gamma3(0.425_dp,3._dp,5._dp,'antikt'))<1e-6_dp) &
    error stop 'Algorithm dependence absent'
  ! Verify conversion to RadISH's non_incl_threeloop convention against
  ! its two hard-coded polynomials, for BOTH quark and gluon representations.
  do k=1,2
    cr=3; if(k==2) cr=4._dp/3
    a1=2*cr; nf=5
    do i=1,2
      r=radii(i)
      if(i==1) then
        legacy=0.0002014631354378078_dp*a1**3 &
          -0.5668546604950793_dp*a1**2*3+299.18580501547285_dp*a1*9 &
          +0.06360185082498425_dp*a1**2*nf/2-161.70214443703713_dp*a1*3*nf/2 &
          +50.82612434499009_dp*a1*(4._dp/3)*nf/2-22.99375663317275_dp*a1*nf**2/4
      else
        legacy=0.050564946123770625_dp*a1**3 &
          -6.313572340177589_dp*a1**2*3+119.6726684453771_dp*a1*9 &
          +1.4719335161834146_dp*a1**2*nf/2-67.05101451465444_dp*a1*3*nf/2 &
          +13.152453101866236_dp*a1*(4._dp/3)*nf/2-7.198933204388625_dp*a1*nf**2/4
      end if
      v=delta_gamma3(r,cr,nf,'antikt')/(128*a1)
      if(abs(v-legacy/(16*a1))>1e-6_dp) error stop 'RadISH normalisation mismatch'
    end do
  end do
end program
