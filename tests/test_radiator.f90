program test_radiator
  use types
  use consts_dp
  use qcd
  use rad_tools
  use radiator_n3ll
  implicit none
  type(process_and_parameters) :: cs
  real(dp) :: aa(4),bbeta(0:3),bbvec(3),c(3),a_small,L,lambda,full(1),diff(2),ratio
  integer :: i,j,k
  call qcd_SetNf(5)
  ! Nonzero independent coefficients exercise every term, not only one process.
  aa=[6._dp,20._dp,130._dp,900._dp]; bbvec=[-7._dp,15._dp,210._dp]
  A=aa(1:3); B=bbvec(1:2)
  bbeta=[beta0,beta1,beta2,0.1_dp]
  L=2.3_dp
  do j=-1,1
    do k=-1,1
      cs%ln_Q2_M2=j*log(4._dp)
      cs%ln_Q2_muR2=k*log(4._dp)
      c=radiator_expansion(L,aa,bbvec,bbeta,cs%ln_Q2_M2,cs%ln_Q2_muR2)
      do i=1,2
        a_small=0.001_dp/2**(i-1)
        cs%as2pi=a_small; cs%alphas_muR=twopi*a_small
        lambda=cs%alphas_muR*beta0*L
        full=Rad([L],cs,order_NNLL)+4*a_small**2 &
          *g4_n3ll(lambda,aa,bbvec,bbeta,cs%ln_Q2_M2,cs%ln_Q2_muR2)
        diff(i)=abs(full(1)-sum(c*[a_small,a_small**2,a_small**3]))
      end do
      ratio=diff(1)/diff(2)
      write(*,'(a,2i3,a,f12.6)') 'scale indices',j,k,' remainder ratio ',ratio
      ! Removing the cubic expansion must leave a fourth-order remainder.
      if(abs(ratio-16._dp)>1._dp) error stop 'Radiator expansion fails order-a^3 test'
    end do
  end do
  if(g4_n3ll(0._dp,aa,bbvec,bbeta,0._dp,0._dp)/=0) error stop 'g4 endpoint is not zero'
  c=radiator_expansion(0._dp,aa,bbvec,bbeta,0._dp,0._dp)
  if(any(c/=0)) error stop 'Expansion endpoint is not zero'
end program
