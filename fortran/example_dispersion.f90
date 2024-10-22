program example_absorption

! example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: profile
   implicit none
   real(dp) :: mHT
   real(dp) :: nu0 = 112265.5949_dp ! Unperturbed line position in cm-1.
   real(dp) :: GammaD = 35.1e-3_dp ! Doppler broadening in cm-1.
   real(dp) :: Gamma0 = 11.3e-3_dp ! Speed-averaged line-width in cm-1.
   real(dp) :: Gamma2 = 0.374e-3_dp ! Speed dependence of the line-width in cm-1.
   real(dp) :: Delta0 = -26.4e-3_dp ! Speed-averaged line-shift in cm-1.
   real(dp) :: Delta2 = 17.8e-3_dp ! Speed dependence of the line-shift in cm-1.
   real(dp) :: NuOptRe = 72.1e-3_dp ! Real part of the Dicke parameter in cm-1.
   real(dp) :: NuOptIm = -16.1e-3_dp ! Imaginary part of the Dicke parameter in cm-1.
   real(dp) :: nu ! Current wavenumber of the computation in cm-1.
   logical :: calculate_dispersion = .true. ! The mHT function will return the dispersion profile.
   nu  = nu0 + 1.0_dp
   mHT = profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,       &
      NuOptIm,nu,calculate_dispersion_opt=calculate_dispersion)
   
   write(*,'(A, F22.15)') "The output of the mHT function (dispersion,"&
      // " Ar-perturbed S(1) 3-0 line in H2):", mHT
end program example_absorption
