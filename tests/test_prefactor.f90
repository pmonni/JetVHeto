program test_prefactor
  use types
  use consts_dp
  use luminosity_n3ll, only: n3ll_prefactor
  implicit none
  real(dp) :: a,b0,b1,L,l0,l1,l2,h1,h2,r1,r2,c(0:3),v,expected,diff(2)
  integer :: i
  l0=7; l1=11; l2=23; h1=3; h2=5; b0=0.61_dp; b1=0.245_dp; L=2
  a=0.1_dp
  v=n3ll_prefactor(l0,l1,l2,a,0._dp,b0,b1,h1,h2)
  expected=l0+a*(l1+h1*l0)+a*a*(l2+h1*l1+h2*l0)
  if(abs(v-expected)>1e-13_dp) error stop 'Spurious higher-order constant in prefactor'
  r1=4*pi*b0*L; r2=8*pi*pi*b1*L
  c(0)=l0
  c(1)=l1+h1*l0
  c(2)=r1*l1+l2+h1*l1+h2*l0
  c(3)=(r1*r1+r2)*l1+2*r1*l2+r1*h1*l1
  do i=1,2
    a=0.001_dp/2**(i-1)
    v=n3ll_prefactor(l0,l1,l2,a,twopi*a*b0*L,b0,b1,h1,h2)
    diff(i)=abs(v-sum(c*[1._dp,a,a*a,a**3]))
  end do
  write(*,*) 'Prefactor fourth-order remainder ratio:',diff(1)/diff(2)
  if(abs(diff(1)/diff(2)-16._dp)>1._dp) error stop 'Prefactor running-coupling order mismatch'
end program
