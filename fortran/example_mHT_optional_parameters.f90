program example_mHT_optional_parameters
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: mHTprofile
   implicit none
   real(dp) :: absorption, dispersion
! Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference: 10.1103/PhysRevA.101.052705)
   real(dp) :: nu0 = 112265.5949_dp ! Unperturbed line position in cm-1.
   real(dp) :: GammaD = 35.1e-3_dp ! Doppler broadening in cm-1.
   real(dp) :: Gamma0 = 11.7e-3_dp ! Speed-averaged line-width in cm-1.
   real(dp) :: Gamma2 = 5.4e-3_dp ! Speed dependence of the line-width in cm-1.
   real(dp) :: Delta0 = 30.5e-3_dp ! Speed-averaged line-shift in cm-1.
   real(dp) :: Delta2 = 12.4e-3_dp ! Speed dependence of the line-shift in cm-1.
   real(dp) :: NuOptRe = 38.0e-3_dp ! Real part of the Dicke parameter in cm-1.
   real(dp) :: NuOptIm = -17.5e-3_dp ! Imaginary part of the Dicke parameter in cm-1.
   real(dp) :: nu ! Current wavenumber of the computation in cm-1.
   real(dp) :: Ylm ! Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
   real(dp) :: Xlm ! Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
   real(dp) :: alpha ! perturber-to-absorber mass ratio, dimensionless.
   logical  :: calculate_dispersion  ! The mHT function will return either the dispersion or absorption profile.
   
   nu  = nu0 + 1.0_dp
   absorption = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,     &
      NuOptRe,NuOptIm,nu)
   calculate_dispersion = .true.
   dispersion = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,     &
      NuOptRe,NuOptIm,nu,calculate_dispersion_opt=calculate_dispersion)
   write(*, '(A, 25X, 2F22.15)') "The output of the mHT function "     &
      // "(He-perturbed S(1) 3-0 line in H2):", absorption, dispersion
! optional parameters
   Ylm = 1.0e-3_dp
   Xlm = 0.5e-3_dp
   alpha = 2.0_dp

   absorption = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,     &
      NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha)
   dispersion = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,     &
      NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)

   write(*, '(A, 2F22.15)') "The output of the mHT function with "     &
      // "optional parameters (He-perturbed S(1) 3-0 line in H2):",    &
      absorption, dispersion
end program example_mHT_optional_parameters
