      program example_mHT_optional_parameters
c Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference: 10.1103/PhysRevA.101.052705)
      implicit none
      double precision :: absorption,dispersion, mHTprofile
      double precision :: nu0 = 112265.5949d0 ! Unperturbed line position in cm-1.
      double precision :: GammaD = 35.1d-3 ! Doppler broadening in cm-1.
      double precision :: Gamma0 = 11.7d-3 ! Speed-averaged line-width in cm-1.
      double precision :: Gamma2 = 5.4d-3 ! Speed dependence of the line-width in cm-1.
      double precision :: Delta0 = 30.5d-3 ! Speed-averaged line-shift in cm-1.
      double precision :: Delta2 = 12.4d-3 ! Speed dependence of the line-shift in cm-1.
      double precision :: NuOptRe = 38.0d-3 ! Real part of the Dicke parameter in cm-1.
      double precision :: NuOptIm = -17.5d-3 ! Imaginary part of the Dicke parameter in cm-1.
      double precision :: nu ! Current wavenumber of the computation in cm-1.
      double precision :: Ylm = 0d0 ! Imaginary part of the 1st order (Rosenkranz) line mixing coefficients
      double precision :: Xlm = 0d0 ! Real part of the 1st order (Rosenkranz) line mixing coefficients
      double precision :: alpha = 10d0 ! perturber-to-absorber mass ratio
      logical :: calculate_dispersion=.false. ! Only the absorption profile will be calculated
   
      nu  = nu0 + 1d0
      absorption = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,
     &            NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)
      calculate_dispersion = .true.
      dispersion = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,
     &            NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)
      write(*, 101) absorption, dispersion
 101  format('The output of the mHT function (He-perturbed S(1)',
     &       ' 3-0 line in H2):',25X, 2F22.15)
c optional parameters
      Ylm = 1.0d-3
      Xlm = 0.5d-3
      alpha = 2.0d0
      calculate_dispersion = .false.
      
      absorption = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,
     &            NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)
     
      calculate_dispersion = .true.
      dispersion = mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,
     &            NuOptRe,NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)

      write(*, 102) absorption, dispersion
 102  format('The output of the mHT function with optional parameters ',
     &       '(He-perturbed S(1) 3-0 line in H2):', 2F22.15)

      end program example_mHT_optional_parameters
