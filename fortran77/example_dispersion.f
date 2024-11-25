      program example_dispersion
c Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference: 10.1103/PhysRevA.101.052705)
      implicit none
      real*8 :: mHT, profile
      real*8 :: nu0 = 112265.5949d0 ! Unperturbed line position in cm-1.
      real*8 :: GammaD = 35.1d-3 ! Doppler broadening in cm-1.
      real*8 :: Gamma0 = 11.7d-3 ! Speed-averaged line-width in cm-1.
      real*8 :: Gamma2 = 5.4d-3 ! Speed dependence of the line-width in cm-1.
      real*8 :: Delta0 = 30.5d-3 ! Speed-averaged line-shift in cm-1.
      real*8 :: Delta2 = 12.4d-3 ! Speed dependence of the line-shift in cm-1.
      real*8 :: NuOptRe = 38.0d-3 ! Real part of the Dicke parameter in cm-1.
      real*8 :: NuOptIm = -17.7d-3 ! Imaginary part of the Dicke parameter in cm-1.
      real*8 :: nu ! Current wavenumber of the computation in cm-1.
      real*8 :: Ylm = 0d0 ! Imaginary part of the 1st order (Rosenkranz) line mixing coefficients
      real*8 :: Xlm = 0d0 ! Real part of the 1st order (Rosenkranz) line mixing coefficients
      real*8 :: alpha = 10d0 ! perturber-to-absorber mass ratio
      logical :: calculate_dispersion = .true. ! The mHT function will return the dispersion profile.

      nu  = nu0 + 1.0d0
      mHT = profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,
     &      NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)
   
      write(*,101) mHT
 101  format('The output of the mHT function (dispersion, He-perturbed',
     &   ' S(1) 3-0 line in H2):', F22.15)
      end program example_dispersion