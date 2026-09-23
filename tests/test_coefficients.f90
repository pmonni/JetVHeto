program test_coefficients
  use hoppet_v1
  use coefficient_functions_n3ll, only: InitCoeffMatrix
  use coefficient_functions_VH, only: InitOldCoeffMatrix => InitCoeffMatrix
  use coefficient_functions_deltaptj
  use convolution_communicator
  use, intrinsic :: ieee_arithmetic
  implicit none
  type(grid_def) :: grid
  type(split_mat) :: c1,c2,g1,delta,old_c1
  real(dp), pointer :: pdf(:,:),result(:,:)
  real(dp) :: x,norm,previous,left,right,jump,pieces(4)
  character(len=5), parameter :: channels(6)=[character(len=5)::'qq','qqbar','qqp','gg','gq','qg']
  integer, parameter :: selectors(4)=[cc_REAL,cc_VIRT,cc_REALVIRT,cc_DELTA]
  integer :: i,j
  character(len=2048) :: data_directory
  call get_command_argument(1,data_directory)
  if(len_trim(data_directory)==0) error stop 'Supply the TMDs_ptj data directory'
  call qcd_SetNf(5)
  call InitGridDef(grid,0.25_dp,8._dp,order=-5)
  call InitCoeffMatrix(grid,c1,c2,g1)
  call InitOldCoeffMatrix(grid,old_c1)
  call AllocPDF(grid,pdf)
  call AllocPDF(grid,result)
  pdf=0
  do i=0,grid%ny
    x=exp(-grid%dy*i)
    pdf(i,0)=x**0.2_dp*(1-x)**5
    do j=1,5
      pdf(i,j)=x**0.3_dp*(1-x)**3/j
      pdf(i,-j)=0.2_dp*pdf(i,j)
    end do
  end do
  result=c1 .conv. pdf
  if(.not.all(ieee_is_finite(result))) error stop 'Nonfinite C1 convolution'
  result=result-(old_c1 .conv. pdf)
  if(maxval(abs(result))>1e-12_dp) error stop 'C1 differs from original JetVHeto'
  result=c2 .conv. pdf
  if(.not.all(ieee_is_finite(result))) error stop 'Nonfinite C2 convolution'
  write(*,*) 'TMD C2 convolution norm:',sum(abs(result))
  previous=0
  do j=1,2
    call InitCoeffMatrixDeltaPtj(grid,delta,0.2_dp*j,trim(data_directory))
    result=delta .conv. pdf
    if(.not.all(ieee_is_finite(result))) error stop 'Nonfinite delta-C2 convolution'
    norm=sum(abs(result))
    if(norm<=0 .or. norm==previous) error stop 'Absent radius dependence'
    write(*,*) 'delta-C2 convolution radius,norm:',0.2_dp*j,norm
    previous=norm
    call Delete(delta)
  end do
  do j=1,6
    do i=1,4
      pieces(i)=evaluate_delta_c2(trim(channels(j)),0.6_dp,selectors(i),0.4_dp,trim(data_directory))
    end do
    if(.not.all(ieee_is_finite(pieces))) error stop 'Nonfinite distribution component'
    if(abs(pieces(1)+pieces(2)-pieces(3))>1e-9_dp*max(1._dp,abs(pieces(3)))) &
      error stop 'Real/virtual decomposition mismatch'
    left=evaluate_delta_c2(trim(channels(j)),0.8_dp-1e-9_dp,cc_REALVIRT,0.4_dp,trim(data_directory))
    right=evaluate_delta_c2(trim(channels(j)),0.8_dp+1e-9_dp,cc_REALVIRT,0.4_dp,trim(data_directory))
    jump=abs(right-left)/max(1._dp,abs(left),abs(right))
    write(*,*) 'Threshold jump diagnostic ',trim(channels(j)),jump
    if(.not.ieee_is_finite(jump)) error stop 'Nonfinite asymptotic threshold'
  end do
  call Delete(old_c1)
  call Delete(c1); call Delete(c2); call Delete(g1)
  call Delete(pdf); call Delete(result); call Delete(grid)
end program
