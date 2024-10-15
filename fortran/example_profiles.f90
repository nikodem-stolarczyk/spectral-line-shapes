program example_profiles
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: profile
   implicit none
   real(dp), allocatable :: nu_tabulated(:) ! Frequency grid in cm-1.
   real(dp), allocatable :: absorption_tabulated(:), dispersion_tabulated(:)
   real(dp) :: nu0 ! Unperturbed line position in cm-1.
   real(dp) :: GammaD ! Doppler broadening in cm-1.
   real(dp) :: Gamma0_He, Gamma0_Ar ! Speed-averaged line-width in cm-1.
   real(dp) :: Gamma2_He, Gamma2_Ar ! Speed dependence of the line-width in cm-1.
   real(dp) :: Delta0_He, Delta0_Ar ! Speed-averaged line-shift in cm-1.
   real(dp) :: Delta2_He, Delta2_Ar ! Speed dependence of the line-shift in cm-1.
   real(dp) :: NuOptRe_He, NuOptRe_Ar ! Real part of the Dicke parameter in cm-1.
   real(dp) :: NuOptIm_He, NuOptIm_Ar ! Imaginary part of the Dicke parameter in cm-1.
   real(dp) :: Ylm ! Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
   real(dp) :: Xlm ! Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
   real(dp) :: alpha_He, alpha_Ar ! perturber-to-absorber mass ratio, dimensionless.
   logical  :: calculate_dispersion ! The mHT function will return either the dispersion or absorption profile.

   integer(int32) :: unit_ = 10 ! File unit number (fixed at 10)
   integer(int32) :: point ! A running index for iteration
   integer(int32) :: number_of_points ! Number of points on the frequency grid
   real(dp) :: nu_step ! Frequency step (in cm^(-1))
   
! example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
   nu0        = 112265.5949_dp
   GammaD     = 35.1e-3_dp
   Gamma0_Ar  = 11.3e-3_dp
   Gamma2_Ar  = 0.374e-3_dp
   Delta0_Ar  =-26.4e-3_dp
   Delta2_Ar  = 17.8e-3_dp
   NuOptRe_Ar = 72.1e-3_dp
   NuOptIm_Ar =-16.1e-3_dp
! optional parameters
   Ylm      = 1.0e-3_dp
   Xlm      = 0.5e-3_dp
   alpha_Ar = 20.0_dp
   calculate_dispersion = .true.

! Generate the mHT profile from -5 GammaD to +5 GammaD and save to an external file
   write(*,*) "Generating the mHT profile for the Ar-perturbed 3-0 S(1) line in H2..."
   number_of_points = 1001
   
   allocate(nu_tabulated(number_of_points))
   allocate(absorption_tabulated(number_of_points))
   allocate(dispersion_tabulated(number_of_points))

   nu_step = (10*GammaD) / real(number_of_points - 1, kind = dp)

   open(unit=unit_, file='mHT_profile_H2Ar.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency         Real(mHT)             Imag(mHT)'
   
   do point = 1, 1001
      nu_tabulated(point) = nu0 - 5*GammaD + (point-1) * nu_step
      absorption_tabulated(point) = profile(nu0,GammaD,Gamma0_Ar,      &
         Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,          &
         nu_tabulated(point),Ylm,Xlm,alpha_Ar)
      dispersion_tabulated(point) = profile(nu0,GammaD,Gamma0_Ar,      &
         Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,          &
         nu_tabulated(point),Ylm,Xlm,alpha_Ar,calculate_dispersion_opt=&
         calculate_dispersion)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         absorption_tabulated(point), dispersion_tabulated(point)

   end do

   close(unit_)
   write(*,*) "The result has been saved to mHT_profile_H2Ar.txt"
   
! example parameters of the S(1) 3-0 line of H2 perturbed by He (reference 10.1103/PhysRevA.101.052705)
   write(*,*) "Generating the mHT profile for the He-perturbed 3-0 S(1) line in H2..."
   Gamma0_He    = 11.7e-3_dp
   Gamma2_He    = 5.4e-3_dp
   Delta0_He  = 30.5e-3_dp
   Delta2_He  = 12.4e-3_dp
   NuOptRe_He = 38.0e-3_dp
   NuOptIm_He =-17.5e-3_dp
   alpha_He = 2.0_dp

   open(unit=unit_, file='mHT_profile_H2He.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency         Real(mHT)             Imag(mHT)'

   do point = 1, 1001
      absorption_tabulated(point) = profile(nu0,GammaD,Gamma0_He,      &
         Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,          &
         nu_tabulated(point),Ylm,Xlm,alpha_He)
      dispersion_tabulated(point) = profile(nu0,GammaD,Gamma0_He,      &
         Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,          &
         nu_tabulated(point),Ylm,Xlm,alpha_He,calculate_dispersion)
      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         absorption_tabulated(point), dispersion_tabulated(point)
   end do

   close(unit_)
   write(*,*) "The result has been saved to mHT_profile_H2He.txt"

end program example_profiles
