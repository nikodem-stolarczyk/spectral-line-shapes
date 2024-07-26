program example_mHT_optional_parameters
!----------------------------------------------------------------------!
! This simple program provides the absorption and dispersion profile
! calls specific for the S(1) 3-0 line of H2 perturbed by Ar
! (reference 10.1063/5.0139229)
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: profile
!----------------------------------------------------------------------!
   implicit none
!----------------------------------------------------------------------!
   complex(dp) :: mHT
      ! Spectral line shape function (in cm); the real part determines
      ! the absorption profile, while its imaginary part provides the
      ! dispersion
   real(dp) :: nu0 = 112265.5949_dp
      ! Central frequency of the transition (in cm^(-1))
   real(dp) :: GamD = 35.1e-3_dp
      ! The Doppler width (in cm^(-1))
   real(dp) :: Gam0 = 11.3e-3_dp
      ! The speed-averaged pressure broadening coefficient (in cm^(-1))
   real(dp) :: Gam2 = 0.374e-3_dp
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: Shift0 = -26.4e-3_dp
      ! The speed-averaged pressure shift coefficient (in cm^(-1))
   real(dp) :: Shift2 = 17.8e-3_dp
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: NuOptRe = 72.1e-3_dp
      ! Real part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: NuOptIm = -16.1e-3_dp
      ! Imaginary part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: nu
      ! Frequency (in cm^(-1))
   real(dp) :: Ylm
      ! Imaginary part of the line mixing coefficient (dimensionless)
   real(dp) :: Xlm
      ! Real part of the line mixing coefficient (dimensionless)
   real(dp) :: alpha
      ! perturber-to-absorber mass ratio (dimensionless)
!----------------------------------------------------------------------!
   nu  = nu0 + 1.0_dp
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu)
   
   write(*, '(A, 25X, 2F22.15)') "The output of the mHT function "     &
      // "(Ar-perturbed S(1) 3-0 line in H2):", real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
! Add the optional parameters
!----------------------------------------------------------------------!
   Ylm = 1.0e-3_dp
   Xlm = 0.5e-3_dp
   alpha = 20.0_dp

   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,  &
      Ylm,Xlm,alpha)
   
   write(*, '(A, 2F22.15)') "The output of the mHT function with "     &
      // "optional parameters (Ar-perturbed S(1) 3-0 line in H2):",    &
      real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
end program example_mHT_optional_parameters
