module coefficient_functions_deltaptj
   use types; use consts_dp; use convolution_communicator; use convolution
   use dglap_objects
   use ew_parameters
   use qcd
   use rad_tools
   use interpolation, only: interpolate
   use reader
   use, intrinsic :: ieee_arithmetic
   implicit none
   private

   public :: InitCoeffMatrixDeltaPtj
   public :: evaluate_delta_c2
   character(len=:), allocatable :: boundary_directory
   real(dp) :: jet_radius
   real(dp), allocatable :: x_and_BCQQBARNS_grid(:,:)
   real(dp), allocatable :: x_and_BCQQBARS_grid(:,:)
   real(dp), allocatable :: x_and_BCQQCACF_grid(:,:)
   real(dp), allocatable :: x_and_BCQQCFCF_grid(:,:)
   real(dp), allocatable :: x_and_BCQQCFTF_grid(:,:)
   real(dp), allocatable :: x_and_BCGGCACA_grid(:,:)
   real(dp), allocatable :: x_and_BCGGCATF_grid(:,:)
   real(dp), allocatable :: x_and_BCGGCFTF_grid(:,:)
   real(dp), allocatable :: x_and_BCGQCFCA_grid(:,:)
   real(dp), allocatable :: x_and_BCGQCFCF_grid(:,:)
   real(dp), allocatable :: x_and_BCQGCATF_grid(:,:)
   real(dp), allocatable :: x_and_BCQGCFTF_grid(:,:)

contains

   ! Diagnostic entry point for channel/distribution and endpoint tests.
   ! Like Hoppet's convolution communicator, this module is not thread-safe.
   function evaluate_delta_c2(channel,x,piece,radius,data_directory) result(value)
      character(len=*), intent(in) :: channel,data_directory
      real(dp), intent(in) :: x,radius
      integer, intent(in) :: piece
      real(dp) :: value
      integer :: saved_piece
      if(.not.ieee_is_finite(x) .or. x<=0 .or. x>=1) error stop 'C2: require 0 < x < 1'
      if(.not.ieee_is_finite(radius) .or. radius<0.05_dp .or. radius>1) error stop 'C2: invalid radius'
      if(piece/=cc_REAL .and. piece/=cc_VIRT .and. piece/=cc_REALVIRT .and. &
        piece/=cc_DELTA) error stop 'C2: unknown distribution piece'
      call init_BC(data_directory)
      jet_radius=radius
      saved_piece=cc_piece; cc_piece=piece
      select case(trim(channel))
      case('qq'); value=deltaC2qq(-log(x))
      case('qqbar'); value=deltaC2qqbar(-log(x))
      case('qqp'); value=deltaC2qqp(-log(x))
      case('gg'); value=deltaC2gg(-log(x))
      case('gq'); value=deltaC2gq(-log(x))
      case('qg'); value=deltaC2qg(-log(x))
      case default; error stop 'C2: unknown channel'
      end select
      cc_piece=saved_piece
   end function

   subroutine init_BC(data_directory)
      character(len=*), intent(in) :: data_directory
      if (allocated(boundary_directory)) then
         if(boundary_directory/=trim(data_directory)) error stop 'C2: data directory changed after initialization'
         return
      end if
      boundary_directory=trim(data_directory)
      call load_boundary("BoundaryConditionQQBARNSFinal.txt", x_and_BCQQBARNS_grid)
      call load_boundary("BoundaryConditionQQBARSFinal.txt", x_and_BCQQBARS_grid)
      call load_boundary("BoundaryConditionQQCACFFinal.txt", x_and_BCQQCACF_grid)
      call load_boundary("BoundaryConditionQQCFCFFinal.txt", x_and_BCQQCFCF_grid)
      call load_boundary("BoundaryConditionQQCFTFFinal.txt", x_and_BCQQCFTF_grid)
      call load_boundary("BoundaryConditionGGCACAFinal.txt", x_and_BCGGCACA_grid)
      call load_boundary("BoundaryConditionGGCATFFinal.txt", x_and_BCGGCATF_grid)
      call load_boundary("BoundaryConditionGGCFTFFinal.txt", x_and_BCGGCFTF_grid)
      call load_boundary("BoundaryConditionGQCFCAFinal.txt", x_and_BCGQCFCA_grid)
      call load_boundary("BoundaryConditionGQCFCFFinal.txt", x_and_BCGQCFCF_grid)
      call load_boundary("BoundaryConditionQGCATFFinal.txt", x_and_BCQGCATF_grid)
      call load_boundary("BoundaryConditionQGCFTFFinal.txt", x_and_BCQGCFTF_grid)
   end subroutine init_BC

   !----------------------------------------------------------------------
   function deltaC2qq(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDQQDelta, TMDQQPlusnoBC, TMDQQPlusnoBCAsym, TMDQQPlus_endpoint
      real(dp) :: TMDQQPlusBCQQBARS, TMDQQPlusBCQQBARNS, TMDQQPlusBCQQCACF, TMDQQPlusBCQQCFCF, TMDQQPlusBCQQCFTF
      real(dp) :: BoundaryConditionQQBARNS, BoundaryConditionQQBARS, BoundaryConditionQQCACF, BoundaryConditionQQCFCF, BoundaryConditionQQCFTF
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDQQ_Delta.txt"
      include "TMDQQ_PlusnoBC.txt"
      include "TMDQQ_PlusnoBCAsym.txt"
      include "TMDQQPlusBCQQBARS.txt"
      include "TMDQQPlusBCQQBARNS.txt"
      include "TMDQQPlusBCQQCACF.txt"
      include "TMDQQPlusBCQQCFCF.txt"
      include "TMDQQPlusBCQQCFTF.txt"
      include "TMDQQPlus_endpoint.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCQQBARNS_grid(1,:), x_and_BCQQBARNS_grid(2,:), x, BoundaryConditionQQBARNS, dy_error, interp_order)
            call interpolate(x_and_BCQQBARS_grid(1,:), x_and_BCQQBARS_grid(2,:), x, BoundaryConditionQQBARS, dy_error, interp_order)
            call interpolate(x_and_BCQQCACF_grid(1,:), x_and_BCQQCACF_grid(2,:), x, BoundaryConditionQQCACF, dy_error, interp_order)
            call interpolate(x_and_BCQQCFCF_grid(1,:), x_and_BCQQCFCF_grid(2,:), x, BoundaryConditionQQCFCF, dy_error, interp_order)
            call interpolate(x_and_BCQQCFTF_grid(1,:), x_and_BCQQCFTF_grid(2,:), x, BoundaryConditionQQCFTF, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDQQPlusnoBC
            else
               res = TMDQQPlusnoBCAsym
            end if
            res = res + TMDQQPlusBCQQBARNS*BoundaryConditionQQBARNS
            res = res + TMDQQPlusBCQQBARS*BoundaryConditionQQBARS
            res = res + TMDQQPlusBCQQCACF*BoundaryConditionQQCACF
            res = res + TMDQQPlusBCQQCFCF*BoundaryConditionQQCFCF
            res = res + TMDQQPlusBCQQCFTF*BoundaryConditionQQCFTF
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
         if (one - x > cutoff) then
            ! subtract the singular part of the 1/(1-x)_+ distribution
            res = res -  TMDQQPlus_endpoint/(one-x)
         else
            res = res
         end if
      case(cc_DELTA)
         res = TMDQQDelta
      end select
      ! write(*,*) x, TMDQQPlusnoBC/(one -x), TMDQQPlusBCQQBARNS/(one -x), TMDQQPlusBCQQBARS/(one -x), TMDQQPlusBCQQCACF/(one -x), TMDQQPlusBCQQCFCF/(one -x), TMDQQPlusBCQQCFTF/(one -x)
      ! write(*,*) BoundaryConditionQQBARNS, BoundaryConditionQQBARS, BoundaryConditionQQCACF, BoundaryConditionQQCFCF, BoundaryConditionQQCFTF
      ! write(*,*) x, res

      if (cc_piece /= cc_DELTA) res = res * x
   end function deltaC2qq

  !----------------------------------------------------------------------
  function deltaC2qqbar(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDQQBPlusnoBC, TMDQQBPlusnoBCAsym
      real(dp) :: TMDQQBPlusBCQQBARS, TMDQQBPlusBCQQBARNS
      real(dp) :: BoundaryConditionQQBARNS, BoundaryConditionQQBARS
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDQQB_PlusnoBC.txt"
      include "TMDQQB_PlusnoBCAsym.txt"
      include "TMDQQBPlusBCQQBARS.txt"
      include "TMDQQBPlusBCQQBARNS.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCQQBARNS_grid(1,:), x_and_BCQQBARNS_grid(2,:), x, BoundaryConditionQQBARNS, dy_error, interp_order)
            call interpolate(x_and_BCQQBARS_grid(1,:), x_and_BCQQBARS_grid(2,:), x, BoundaryConditionQQBARS, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDQQBPlusnoBC
            else
               res = TMDQQBPlusnoBCAsym
            end if
            res = res + TMDQQBPlusBCQQBARNS*BoundaryConditionQQBARNS
            res = res + TMDQQBPlusBCQQBARS*BoundaryConditionQQBARS
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
      case(cc_DELTA)
         res = zero
      end select

      if (cc_piece /= cc_DELTA) res = res * x
  end function deltaC2qqbar

  !----------------------------------------------------------------------
  function deltaC2qqp(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDQQPPlusnoBC, TMDQQPPlusnoBCAsym
      real(dp) :: TMDQQPPlusBCQQBARS
      real(dp) :: BoundaryConditionQQBARS
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDQQP_PlusnoBC.txt"
      include "TMDQQP_PlusnoBCAsym.txt"
      include "TMDQQPPlusBCQQBARS.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCQQBARS_grid(1,:), x_and_BCQQBARS_grid(2,:), x, BoundaryConditionQQBARS, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDQQPPlusnoBC
            else
               res = TMDQQPPlusnoBCAsym
            end if
            res = res + TMDQQPPlusBCQQBARS*BoundaryConditionQQBARS
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
      case(cc_DELTA)
         res = zero
      end select

      if (cc_piece /= cc_DELTA) res = res * x
  end function deltaC2qqp

  !----------------------------------------------------------------------
  function deltaC2gq(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDGQPlusnoBC, TMDGQPlusnoBCAsym
      real(dp) :: TMDGQPlusBCGQCFCA, TMDGQPlusBCGQCFCF
      real(dp) :: BoundaryConditionGQCFCA, BoundaryConditionGQCFCF
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDGQ_PlusnoBC.txt"
      include "TMDGQ_PlusnoBCAsym.txt"
      include "TMDGQPlusBCGQCFCA.txt"
      include "TMDGQPlusBCGQCFCF.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCGQCFCF_grid(1,:), x_and_BCGQCFCF_grid(2,:), x, BoundaryConditionGQCFCF, dy_error, interp_order)
            call interpolate(x_and_BCGQCFCA_grid(1,:), x_and_BCGQCFCA_grid(2,:), x, BoundaryConditionGQCFCA, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDGQPlusnoBC
            else
               res = TMDGQPlusnoBCAsym
            end if
            res = res + TMDGQPlusBCGQCFCF*BoundaryConditionGQCFCF
            res = res + TMDGQPlusBCGQCFCA*BoundaryConditionGQCFCA
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
      case(cc_DELTA)
         res = zero
      end select

      if (cc_piece /= cc_DELTA) res = res * x
  end function deltaC2gq

  !----------------------------------------------------------------------
  function deltaC2gg(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDGGPlusnoBC, TMDGGPlusnoBCAsym, TMDGGDelta, TMDGGPlus_endpoint
      real(dp) :: TMDGGPlusBCGGCACA, TMDGGPlusBCGGCATF, TMDGGPlusBCGGCFTF
      real(dp) :: BoundaryConditionGGCACA, BoundaryConditionGGCATF, BoundaryConditionGGCFTF
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDGG_Delta.txt"
      include "TMDGG_PlusnoBC.txt"
      include "TMDGG_PlusnoBCAsym.txt"
      include "TMDGGPlusBCGGCACA.txt"
      include "TMDGGPlusBCGGCATF.txt"
      include "TMDGGPlusBCGGCFTF.txt"
      include "TMDGGPlus_endpoint.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCGGCACA_grid(1,:), x_and_BCGGCACA_grid(2,:), x, BoundaryConditionGGCACA, dy_error, interp_order)
            call interpolate(x_and_BCGGCATF_grid(1,:), x_and_BCGGCATF_grid(2,:), x, BoundaryConditionGGCATF, dy_error, interp_order)
            call interpolate(x_and_BCGGCFTF_grid(1,:), x_and_BCGGCFTF_grid(2,:), x, BoundaryConditionGGCFTF, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDGGPlusnoBC
            else
               res = TMDGGPlusnoBCAsym
            end if
            res = res + TMDGGPlusBCGGCACA*BoundaryConditionGGCACA
            res = res + TMDGGPlusBCGGCATF*BoundaryConditionGGCATF
            res = res + TMDGGPlusBCGGCFTF*BoundaryConditionGGCFTF
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
         if (one - x > cutoff) then
            ! subtract the singular part of the 1/(1-x)_+ distribution
            res = res -  TMDGGPlus_endpoint/(one-x)
         else
            res = res
         end if
      case(cc_DELTA)
         res = TMDGGDelta
      end select

      if (cc_piece /= cc_DELTA) res = res * x
  end function deltaC2gg

  !----------------------------------------------------------------------

   !----------------------------------------------------------------------
   function deltaC2qg(y) result(res)
      real(dp), intent(in) :: y
      real(dp)             :: res
      real(dp)             :: x
      complex(dp) :: Li2x, Li2mx
      real(dp), parameter    :: cutoff=1e-11 ! cutoff to regularise the x->1 singularity in the gg coefficient function
      real(dp) :: TMDQGPlusnoBC, TMDQGPlusnoBCAsym
      real(dp) :: TMDQGPlusBCQGCFTF, TMDQGPlusBCQGCATF
      real(dp) :: BoundaryConditionQGCFTF, BoundaryConditionQGCATF
      real(dp) :: R
      real(dp), parameter :: AsymThreshold=0.8_dp
      !------------------------------------------
      ! The following declarations are necessary for hplog
      integer, parameter :: n1=-1
      integer, parameter :: n2= 1
      integer, parameter :: nw= 4
      integer :: interp_order = 3
      real(dp):: dy_error
      real(dp)    :: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2),&
         & Hr4(n1:n2,n1:n2,n1:n2,n1:n2)
      real(dp)    :: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),&
         & Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      complex(dp) :: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), &
         & Hc4(n1:n2,n1:n2,n1:n2,n1:n2)


      x = exp(-y)
      res = zero

      R = jet_radius

      ! evaluate HPL's using hplog (much faster than Chaplin)
      call hplog(x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2x = Hc2(0,1)
      call hplog(-x,nw,Hc1,Hc2,Hc3,Hc4, &
            &     Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      Li2mx = Hc2(0,1)
      ! TODO: do not use your paths
      include "TMDQG_PlusnoBC.txt"
      include "TMDQG_PlusnoBCAsym.txt"
      include "TMDQGPlusBCQGCFTF.txt"
      include "TMDQGPlusBCQGCATF.txt"
      select case(cc_piece)
      case(cc_REAL,cc_REALVIRT)
         if (one - x > cutoff) then
            call interpolate(x_and_BCQGCATF_grid(1,:), x_and_BCQGCATF_grid(2,:), x, BoundaryConditionQGCATF, dy_error, interp_order)
            call interpolate(x_and_BCQGCFTF_grid(1,:), x_and_BCQGCFTF_grid(2,:), x, BoundaryConditionQGCFTF, dy_error, interp_order)
            if (x < AsymThreshold) then
               res = TMDQGPlusnoBC
            else
               res = TMDQGPlusnoBCAsym
            end if
            res = res + TMDQGPlusBCQGCATF*BoundaryConditionQGCATF
            res = res + TMDQGPlusBCQGCFTF*BoundaryConditionQGCFTF
            res = res / (one-x)
         else
            res = zero
         end if
      end select
      select case(cc_piece)
      case(cc_VIRT,cc_REALVIRT)
      case(cc_DELTA)
         res = zero
      end select

      if (cc_piece /= cc_DELTA) res = res * x
   end function deltaC2qg

  !======================================================================
  !! Initialise a coefficient function matrix, repeating exactly
  !! what's done for a LO unpolarised splitting matrix, with the nf
  !! value that is current from the qcd module.
  subroutine InitCoeffMatrixDeltaPtj(grid, deltaC2, jet_radius_in, data_directory)
    character(len=*), intent(in) :: data_directory
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: deltaC2
    real(dp), intent(in) :: jet_radius_in
    type(grid_conv) :: deltaC2qqS_grid, deltaC2qqbar_grid
    integer :: i
    real(dp):: y, tmp1

    write(*,*) 'Initialising deltaptj coefficient function grids ...'

    call init_BC(data_directory)
    if (.not.ieee_is_finite(jet_radius_in) .or. jet_radius_in<0.05_dp .or. jet_radius_in>1._dp) &
      error stop 'C2: supported jet radius is 0.05 <= R <= 1'
    !-- O(as^2) coefficient functions. The deltaC2qq is non-trivial at this order and we need
    !-- to decompose it consistently in the Hoppet's evolution basis (singlet, non-singlet, valence)
    !deltaC2%loops  = 1
    ! fix the number of active flavours
    deltaC2%nf_int = nf_int

   jet_radius = jet_radius_in

    call cobj_InitSplitLinks(deltaC2)

    !-- Interpret deltaC2qq and deltaC2qqbar as non-singlet components, and
    !-- deltaC2qq', deltaC2qqbar' as singlet components. At as^2 deltaC2qq'=deltaC2qqbar' as
    !-- for the singlet splitting functions. Therefore it seems that all
    !-- the symmetries are preserved. This follows Hoppet's implementation
    !-- of the NLO splitting functions

    !-- start with the +- components of the non-singlet
    !-- PqqV +- PqqbarV, remember that deltaC2qq'=deltaC2qqbar' at as^2
    ! plus component
    call InitGridConv(grid, deltaC2%NS_plus,  deltaC2qq)
    call InitGridConv(grid, deltaC2qqbar_grid,  deltaC2qqbar)
    call AddWithCoeff(deltaC2%NS_plus,  deltaC2qqbar_grid, one)
    call InitGridConv(grid, deltaC2qqS_grid,  deltaC2qqp)
    call AddWithCoeff(deltaC2%NS_plus,  deltaC2qqS_grid, -two)
    ! minus component
    call InitGridConv(grid, deltaC2%NS_minus, deltaC2qq)
    call AddWithCoeff(deltaC2%NS_minus, deltaC2qqbar_grid, -one)

    !-- now the valence component
    !-- PNSminus + nf * (PqqS - PqqbarS), remember that deltaC2qq'=deltaC2qqbar' at as^2
    call InitGridConv(deltaC2%NS_V, deltaC2%NS_minus)

    !-- finally, build the singlet component (eq. 2.4 of Moch et al. 0404111)
    !-- Pqq = PNS_plus + nf*(PqqS + PqqbarS), remember that deltaC2qq'=deltaC2qqbar' at as^2
    call InitGridConv(deltaC2%qq, deltaC2%NS_plus)
    call AddWithCoeff(deltaC2%qq, deltaC2qqS_grid, two*nf)

    ! Clean up
    call Delete(deltaC2qqS_grid)
    call Delete(deltaC2qqbar_grid)

    !-- now build the rest of the singlet
    call InitGridConv(grid, deltaC2%gg, deltaC2gg)
    call InitGridConv(grid, deltaC2%gq, deltaC2gq)
    call InitGridConv(grid, deltaC2%qg, deltaC2qg)
    !-- add the factor of 2*nf in Cqg to sum over all species
    call Multiply(deltaC2%qg, two*nf)

   end subroutine InitCoeffMatrixDeltaPtj


   subroutine load_boundary(filename,values)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: values(:,:)
      type(reader_opts) :: opts_grid
      logical :: exists
      character(len=:), allocatable :: path
      path=boundary_directory//'/BoundaryConditions/'//filename
      inquire(file=path,exist=exists)
      if(.not.exists) then
         write(*,'(a)') 'Missing C2 boundary data: '//path
         error stop 'C2: missing boundary data'
      end if
      call read_file_into_array(path,values,opts_grid,max_col1=1e5_dp)
      if(size(values,1)<2 .or. size(values,2)<3) error stop 'C2: invalid boundary grid shape'
      if(.not.all(ieee_is_finite(values))) error stop 'C2: nonfinite boundary data'
      if(any(values(1,2:)<=values(1,:size(values,2)-1))) error stop 'C2: unordered boundary grid'
   end subroutine

end module coefficient_functions_deltaptj
