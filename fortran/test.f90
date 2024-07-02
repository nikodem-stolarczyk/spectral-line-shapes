program test_mHT
!----------------------------------------------------------------------!
! This simple program tests the subroutines computing the modified
! Hartmann-Tran (mHT) profile for a set of default line-shape parameters
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module
!----------------------------------------------------------------------!
   implicit none
!----------------------------------------------------------------------!
   complex(dp) :: mHT
      ! Spectral line shape function (in cm)
   real(dp) :: nu0     = 0.0_dp
      ! Central frequency of the transition (in cm^(-1))
   real(dp) :: GamD    = 1.0_dp
      ! The Doppler width (in cm^(-1))
   real(dp) :: Gam0    = 1.0_dp
      ! The speed-averaged pressure broadening coefficient (in cm^(-1))
   real(dp) :: Gam2    = 0.1_dp
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: Shift0  = 1.0_dp
      ! The speed-averaged pressure shift coefficient (in cm^(-1))
   real(dp) :: Shift2  = 0.1_dp
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: NuOptRe = 0.1_dp
      ! Real part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: NuOptIm = 0.01_dp
      ! Imaginary part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: nu      = 0.0_dp
      ! Frequency (in cm^(-1))
!----------------------------------------------------------------------!
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu)
   print *, "mHT test: ", mHT
!----------------------------------------------------------------------!
end program test_mHT
