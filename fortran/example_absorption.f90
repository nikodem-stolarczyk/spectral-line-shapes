program example_absorption

! example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: profile
   implicit none
   real(dp) :: mHT
      ! Spectral line shape function (in cm); absorption profile is
      ! obtained as the real part of the profile() function:
      ! "real(profile(...))"
   real(dp) :: nu0 = 112265.5949_dp ! Unperturbed line position in cm-1.
   real(dp) :: GamD = 35.1e-3_dp ! Doppler broadening in cm-1.
   real(dp) :: Gam0 = 11.3e-3_dp ! Speed-averaged line-width in cm-1.
   real(dp) :: Gam2 = 0.374e-3_dp ! Speed dependence of the line-width in cm-1.
   real(dp) :: Shift0 = -26.4e-3_dp ! Speed-averaged line-shift in cm-1.
   real(dp) :: Shift2 = 17.8e-3_dp ! Speed dependence of the line-shift in cm-1.
   real(dp) :: NuOptRe = 72.1e-3_dp ! Real part of the Dicke parameter in cm-1.
   real(dp) :: NuOptIm = -16.1e-3_dp ! Imaginary part of the Dicke parameter in cm-1.
   real(dp) :: nu ! Current wavenumber of the computation in cm-1.
   nu  = nu0 + 1.0_dp
   mHT = real(profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu))
   
   write(*,'(A, F22.15)') "The output of the mHT function (absorption,"&
      // " Ar-perturbed S(1) 3-0 line in H2):", mHT
end program example_absorption
