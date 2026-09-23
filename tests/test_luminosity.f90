! PDF/convolution smoke test, NOT a physical N3LL cross-section benchmark.
program test_luminosity
  use types
  use consts_dp
  use qcd, only: beta0
  use rad_tools
  use pdfs_tools
  use luminosity_n3ll
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(process_and_parameters) :: cs
  character(len=2048) :: directory,pdfname
  real(dp) :: value,base(1),a_saved,lambda
  integer :: i
  call get_command_argument(1,directory)
  call get_command_argument(2,pdfname)
  if(len_trim(pdfname)==0) error stop 'Supply a PDF set for the luminosity smoke test'
  call init_pdfs_from_LHAPDF(trim(pdfname),0)
  call set_process_and_parameters(cs,'pp','H',13000._dp,125._dp,125._dp,125._dp, &
    125._dp,5._dp,0.4_dp,'none','ptj')
  call init_proc(cs)
  call init_luminosity_n3ll(cs%jet_radius,trim(directory))
  a_saved=cs%as2pi
  cs%as2pi=0
  base=lumi_NLL([30._dp],cs)
  value=n3ll_lumi(30._dp,0._dp,cs,H(1),0._dp)
  if(abs(value-base(1))>1e-12_dp*abs(base(1))) error stop 'Zero-coupling luminosity mismatch'
  cs%as2pi=a_saved
  do i=-1,1
    cs%ln_Q2_M2=i*log(4._dp)
    cs%ln_Q2_muR2=cs%ln_Q2_M2
    lambda=cs%alphas_muR*beta0*log(125._dp/30)
    ! hard2 is deliberately set to zero to test assembly, not physics.
    value=n3ll_lumi(30._dp,lambda,cs,H(1),0._dp)
    if(.not.ieee_is_finite(value)) error stop 'Nonfinite assembled luminosity'
    write(*,*) 'Luminosity assembly smoke test, scale index:',i,' value:',value
  end do
end program
