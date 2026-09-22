module ew_parameters
  use types; use consts_dp
  implicit none
  private
  
  real(dp), parameter, public :: MZ=91.1876_dp
  ! PFM hack 01/2015
  real(dp), parameter, public :: sinwsq =  0.2226459_dp ! same as MCFM
!  real(dp), parameter, public :: sinwsq =  0.75_dp ! same as MCFM for Z' (MZ'=160.77 GeV)
!  real(dp), parameter, public :: sinwsq =  0.868887504565376_dp ! same as MCFM for Z' (MZ'=222 GeV)
  real(dp), parameter, public :: invGev2_to_nb = 389379.323_dp
  real(dp), parameter, public :: higgs_vev = 246.21846_dp ! from MCFM 
  real(dp), parameter, public :: GF_GeVm2 = 0.116639e-04_dp ! from MCFM
  real(dp), parameter, public :: gv2_ga2_u = (half-four/three*sinwsq)**2 + one/four
  real(dp), parameter, public :: gv2_ga2_d = (-half+two/three*sinwsq)**2 + one/four
  real(dp), parameter, public :: gv2_ga2(6) = (/gv2_ga2_d, gv2_ga2_u,&
                                              & gv2_ga2_d, gv2_ga2_u,&
                                              & gv2_ga2_d, gv2_ga2_u /)
  real(dp), parameter, public :: mb = 4.75   ! from PDG 2012
  real(dp), parameter, public :: mt = 172.5  ! from PDG 2012

!  real(dp), parameter, public :: mb = 4.65   ! from PDG 2012
!  real(dp), parameter, public :: mt = 173.5  ! from PDG 2012
end module ew_parameters

