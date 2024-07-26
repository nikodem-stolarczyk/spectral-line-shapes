program example_profiles
!----------------------------------------------------------------------!
! This simple program provides two spectral line shape profiles
! (He- and Ar-perturbed 3-0 S(1) line in H2) on a frequency grid.
! The output is saved to two text files.
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use spectral_module, only: profile
!----------------------------------------------------------------------!
   implicit none
!----------------------------------------------------------------------!
! mHT test:
!----------------------------------------------------------------------!
   real(dp), allocatable :: nu_tabulated(:)
      ! Frequency grid (in cm^(-1))
   complex(dp), allocatable :: mHT_tabulated(:)
      ! A tabulated spectral line shape function (in cm). The real part
      ! determines the absorption profile, while its imaginary part
      ! provides the dispersion
   real(dp) :: nu0
      ! Central frequency of the transition (in cm^(-1))
   real(dp) :: GamD
      ! The Doppler width (in cm^(-1))
   real(dp) :: Gam0_He, Gam0_Ar
      ! The speed-averaged pressure broadening coefficient (in cm^(-1))
   real(dp) :: Gam2_He, Gam2_Ar
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: Shift0_He, Shift0_Ar
      ! The speed-averaged pressure shift coefficient (in cm^(-1))
   real(dp) :: Shift2_He, Shift2_Ar
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: NuOptRe_He, NuOptRe_Ar
      ! Real part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: NuOptIm_He, NuOptIm_Ar
      ! Imaginary part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: Ylm
      ! Imaginary part of the line mixing coefficient (dimensionless)
   real(dp) :: Xlm
      ! Real part of the line mixing coefficient (dimensionless)
   real(dp) :: alpha_He, alpha_Ar
      ! perturber-to-absorber mass ratio (dimensionless)
!----------------------------------------------------------------------!
   integer(int32) :: unit_ = 10
      ! File unit number (fixed at 10)
   integer(int32) :: point
      ! A running index for iteration
   integer(int32) :: number_of_points
      ! Number of points on the frequency grid
   real(dp) :: nu_step
      ! Frequency step (in cm^(-1))
!----------------------------------------------------------------------!
! mHT test: example parameters of the S(1) 3-0 line of H2
! perturbed by Ar (reference 10.1063/5.0139229)
!----------------------------------------------------------------------!
   nu0        = 112265.5949_dp
   GamD       = 35.1e-3_dp
   Gam0_Ar    = 11.3e-3_dp
   Gam2_Ar    = 0.374e-3_dp
   Shift0_Ar  =-26.4e-3_dp
   Shift2_Ar  = 17.8e-3_dp
   NuOptRe_Ar = 72.1e-3_dp
   NuOptIm_Ar =-16.1e-3_dp
!----------------------------------------------------------------------!
! optional parameters
!----------------------------------------------------------------------!
   Ylm      = 1.0e-3_dp
   Xlm      = 0.5e-3_dp
   alpha_Ar = 20.0_dp
!----------------------------------------------------------------------!
! Generate the mHT profile from -5 GammaD to +5 GammaD
! and save to an external file
!----------------------------------------------------------------------!
   write(*,*) "Generating the mHT profile for the Ar-perturbed 3-0 S(1) line in H2..."
   number_of_points = 1001
   
   allocate(nu_tabulated(number_of_points))
   allocate(mHT_tabulated(number_of_points))

   nu_step = (10*GamD) / real(number_of_points - 1, kind = dp)

   open(unit=unit_, file='mHT_profile_H2Ar.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency         Real(mHT)             Imag(mHT)'
   
   do point = 1, 1001
      nu_tabulated(point) = nu0 - 5*GamD + (point-1) * nu_step
      mHT_tabulated(point) = profile(nu0,GamD,Gam0_Ar,Gam2_Ar,         &
         Shift0_Ar,Shift2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu_tabulated(point),&
         Ylm,Xlm,alpha_Ar)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         real(mHT_tabulated(point)), aimag(mHT_tabulated(point))

   end do

   close(unit_)
   write(*,*) "The result has been saved to mHT_profile_H2Ar.txt"
!----------------------------------------------------------------------!
! mHT test: example parameters of the S(1) 3-0 line of H2
! perturbed by He (reference 10.1103/PhysRevA.101.052705)
!----------------------------------------------------------------------!
   write(*,*) "Generating the mHT profile for the He-perturbed 3-0 S(1) line in H2..."
   Gam0_He    = 11.7e-3_dp
   Gam2_He    = 5.4e-3_dp
   Shift0_He  = 30.5e-3_dp
   Shift2_He  = 12.4e-3_dp
   NuOptRe_He = 38.0e-3_dp
   NuOptIm_He =-17.5e-3_dp

   alpha_He = 2.0_dp

   open(unit=unit_, file='mHT_profile_H2He.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency         Real(mHT)             Imag(mHT)'

   do point = 1, 1001
      mHT_tabulated(point) = profile(nu0,GamD,Gam0_He,Gam2_He,         &
         Shift0_He,Shift2_He,NuOptRe_He,NuOptIm_He,nu_tabulated(point),&
         Ylm,Xlm,alpha_He)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         real(mHT_tabulated(point)), aimag(mHT_tabulated(point))

   end do

   close(unit_)
   write(*,*) "The result has been saved to mHT_profile_H2He.txt"
!----------------------------------------------------------------------!
end program example_profiles
