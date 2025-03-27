      program example_profiles
      implicit none
      integer, parameter :: number_of_points = 1001 ! Number of points on the frequency grid
      double precision :: nu_tabulated(number_of_points) ! Frequency grid in cm-1.
      double precision :: absorption_tabulated(number_of_points),
     &          dispersion_tabulated(number_of_points)
      double precision :: mHTprofile
      double precision :: nu0 ! Unperturbed line position in cm-1.
      double precision :: GammaD ! Doppler broadening in cm-1.
      double precision :: Gamma0_He, Gamma0_Ar ! Speed-averaged line-width in cm-1.
      double precision :: Gamma2_He, Gamma2_Ar ! Speed dependence of the line-width in cm-1.
      double precision :: Delta0_He, Delta0_Ar ! Speed-averaged line-shift in cm-1.
      double precision :: Delta2_He, Delta2_Ar ! Speed dependence of the line-shift in cm-1.
      double precision :: NuOptRe_He, NuOptRe_Ar ! Real part of the Dicke parameter in cm-1.
      double precision :: NuOptIm_He, NuOptIm_Ar ! Imaginary part of the Dicke parameter in cm-1.
      double precision :: Ylm ! Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
      double precision :: Xlm ! Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
      double precision :: alpha_He, alpha_Ar ! perturber-to-absorber mass ratio, dimensionless.
      logical :: calculate_dispersion, not_calculate_dispersion ! The mHT function will return either the dispersion or absorption profile.

      integer :: unit_ = 10 ! File unit number (fixed at 10)
      integer :: point ! A running index for iteration
      double precision  :: nu_step ! Frequency step (in cm^(-1))
   
c example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
      nu0        = 112265.5949d0
      GammaD     = 35.1d-3
      Gamma0_Ar  = 11.3d-3
      Gamma2_Ar  = 0.374d-3
      Delta0_Ar  =-26.4d-3
      Delta2_Ar  = 17.8d-3
      NuOptRe_Ar = 72.1d-3
      NuOptIm_Ar =-16.1d-3
      Ylm        = 1.0d-3
      Xlm        = 0.5d-3
      alpha_Ar   = 20d0
      calculate_dispersion = .true.
      not_calculate_dispersion = .false.

c Generate the mHT profile from -5 GammaD to +5 GammaD and save to an external file
      write(*,101)
 101  format('Generating the mHT profile for the Ar-perturbed',
     &       ' 3-0 S(1) line in H2...')
      nu_step = (10*GammaD) / dfloat(number_of_points - 1)

      open(unit=unit_, file='mHT_profile_H2Ar.txt', status='replace',
     &     action='write')

      write(unit_, '(A11,9X,A9,13X,A9)') '# Frequency', 'Real(mHT)',
     &     'Imag(mHT)'
   
      do point = 1, number_of_points
         nu_tabulated(point) = nu0 - 5*GammaD + (point-1) * nu_step
         absorption_tabulated(point) = mHTprofile(nu0,GammaD,Gamma0_Ar,
     &      Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,
     &      nu_tabulated(point),Ylm,Xlm,alpha_Ar,
     &      not_calculate_dispersion)
         dispersion_tabulated(point) = mHTprofile(nu0,GammaD,Gamma0_Ar,
     &      Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,
     &      nu_tabulated(point),Ylm,Xlm,alpha_Ar,calculate_dispersion)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),
     &   absorption_tabulated(point), dispersion_tabulated(point)

      end do

      close(unit_)
      write(*,102)
 102  format('The result has been saved to mHT_profile_H2Ar.txt')
   
c example parameters of the S(1) 3-0 line of H2 perturbed by He (reference 10.1103/PhysRevA.101.052705)
      write(*,103)
 103  format('Generating the mHT profile for the He-perturbed 3-0 S(1)', 
     &       ' line in H2...')
      Gamma0_He  = 11.7d-3
      Gamma2_He  = 5.4d-3
      Delta0_He  = 30.5d-3
      Delta2_He  = 12.4d-3
      NuOptRe_He = 38.0d-3
      NuOptIm_He =-17.5d-3
      alpha_He = 2d0

      open(unit=unit_, file='mHT_profile_H2He.txt', status='replace',
     &     action='write')

      write(unit_, '(A11,9X,A9,13X,A9)') '# Frequency', 'Real(mHT)',
     &     'Imag(mHT)'

      do point = 1, number_of_points
         absorption_tabulated(point) = mHTprofile(nu0,GammaD,Gamma0_He,
     &      Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,
     &      nu_tabulated(point),Ylm,Xlm,alpha_He,
     &      not_calculate_dispersion)
         dispersion_tabulated(point) = mHTprofile(nu0,GammaD,Gamma0_He,
     &      Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,
     &      nu_tabulated(point),Ylm,Xlm,alpha_He,calculate_dispersion)
      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),
     &   absorption_tabulated(point), dispersion_tabulated(point)
      end do

      close(unit_)
      write(*,104)
 104  format('The result has been saved to mHT_profile_H2He.txt')

      end program example_profiles
