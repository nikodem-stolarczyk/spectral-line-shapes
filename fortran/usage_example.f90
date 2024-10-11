program test_mHT
!----------------------------------------------------------------------!
! This simple program tests the subroutines computing the complex
! probability function (CPF and the modified Hartmann-Tran (mHT) profile
! for a set of default line-shape parameters.
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use cpf_module, only: cpf_accurate, cpf_fast
   use spectral_module, only: profile
!----------------------------------------------------------------------!
   implicit none
!----------------------------------------------------------------------!
! CPF test:
!----------------------------------------------------------------------!
   real(dp)    :: cpf_input_real = 1.0_dp
      ! Example real argument of the input for the CPF
   real(dp)    :: cpf_input_imag = 1.0_dp
      ! Example imaginary argument of the input for the CPF
   complex(dp) :: cpf_accurate_output
      ! Example output of cpf_accurate(cpf_input_real, cpf_input_imag)
   complex(dp) :: cpf_fast_output
      ! Example output of cpf_fast(cpf_input_real, cpf_input_imag)
!----------------------------------------------------------------------!
! mHT test:
!----------------------------------------------------------------------!
   complex(dp) :: mHT
      ! Spectral line shape function (in cm)
   real(dp) :: nu0
      ! Central frequency of the transition (in cm^(-1))
   real(dp) :: GamD
      ! The Doppler width (in cm^(-1))
   real(dp) :: Gam0
      ! The speed-averaged pressure broadening coefficient (in cm^(-1))
   real(dp) :: Gam2
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: Shift0
      ! The speed-averaged pressure shift coefficient (in cm^(-1))
   real(dp) :: Shift2
      ! The quadratic speed-dependence coefficient for pressure broadening
      ! (in cm^(-1))
   real(dp) :: NuOptRe
      ! Real part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: NuOptIm
      ! Imaginary part of the complex Dicke parameter (in cm^(-1))
   real(dp) :: Sw
      ! Statistical weight (dimensionless)
   real(dp) :: Ylm
      ! Imaginary part of the line mixing coefficient (dimensionless)
   real(dp) :: Xlm
      ! Real part of the line mixing coefficient (dimensionless)
   real(dp) :: alpha
      ! Mass ratio 
   real(dp) :: nu
      ! Frequency (in cm^(-1))
!----------------------------------------------------------------------!
   integer(int32) :: unit_ = 10
      ! File unit number (fixed at 10)
   integer(int32) :: point
      ! A running index for iteration
   integer(int32) :: number_of_points
      ! Number of points on the frequency grid
   real(dp) :: nu_step
      ! Frequency step (in cm^(-1))
   real(dp), allocatable :: nu_tabulated(:)
      ! Frequency grid (in cm^(-1))
   complex(dp), allocatable :: mHT_tabulated(:)
      ! A tabulated spectral line shape function (in cm)
!----------------------------------------------------------------------!
! CPF test: cpf_accurate and cpf_fast
!----------------------------------------------------------------------!
   cpf_input_real = 1.0_dp
   cpf_input_imag = 1.0_dp

   cpf_accurate_output = cpf_accurate(cpf_input_real, cpf_input_imag)
   cpf_fast_output = cpf_fast(cpf_input_real, cpf_input_imag)

   write(*, '(A, 2F22.15)') 'The output of the cpf_accurate function: ',&
      real(cpf_accurate_output), aimag(cpf_accurate_output)
   write(*, '(A, 2F22.15)') 'The output of the cpf_fast function:     ',&
      real(cpf_fast_output), aimag(cpf_fast_output)
   write(*,*)
!----------------------------------------------------------------------!
! mHT test: example parameters of the S(1) 3-0 line of H2
! perturbed by Ar (reference 10.1063/5.0139229)
!----------------------------------------------------------------------!
   nu0     = 112265.5949_dp
   GamD    = 35.1e-3_dp
   Gam0    = 11.3e-3_dp
   Gam2    = 0.374e-3_dp
   Shift0  =-26.4e-3_dp
   Shift2  = 17.8e-3_dp
   NuOptRe = 72.1e-3_dp
   NuOptIm =-16.1e-3_dp
   nu      = nu0 + 1.0_dp
   
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu)
   
   write(*,'(A, 2F22.15)') "The output of the mHT function "           &
      // "(Ar-perturbed S(1) 3-0 line in H2):", real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
! mHT test: the same line but with additional (optional) parameters
!----------------------------------------------------------------------!
   Sw    = 1.0_dp
   Ylm   = 1.0e-3_dp
   Xlm   = 0.5e-3_dp
   alpha = 20.0_dp
   
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,  &
      Sw,Ylm,Xlm,alpha)

   write(*,'(A,22X,2F22.15)') "The same function with optional "       &
      // "parameters: ", real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
! Generate the mHT profile from -5 GammaD to +5 GammaD
! and save to an external file
!----------------------------------------------------------------------!
   number_of_points = 1001
   
   allocate(nu_tabulated(number_of_points))
   allocate(mHT_tabulated(number_of_points))

   nu_step = (10*GamD) / real(number_of_points - 1, kind = dp)

   open(unit=unit_, file='mHT_profile_H2Ar.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency       Real(mHT)             Imag(mHT)'
   
   do point = 1, 1001
      nu_tabulated(point) = nu0 - 5*GamD + (point-1) * nu_step
      mHT_tabulated(point) = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2, &
         NuOptRe,NuOptIm,nu_tabulated(point),Sw,Ylm,Xlm,alpha)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         real(mHT_tabulated(point)), aimag(mHT_tabulated(point))

   end do

   close(unit_)
!----------------------------------------------------------------------!
! mHT test: example parameters of the S(1) 3-0 line of H2
! perturbed by He (reference 10.1103/PhysRevA.101.052705)
!----------------------------------------------------------------------!
   Gam0    = 11.7e-3_dp
   Gam2    = 5.4e-3_dp
   Shift0  = 30.5e-3_dp
   Shift2  = 12.4e-3_dp
   NuOptRe = 38.0e-3_dp
   NuOptIm =-17.5e-3_dp
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu)
   write(*, '(A, 2F22.15)') "The output of the mHT function "          &
      // "(He-perturbed S(1) 3-0 line in H2):", real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
! mHT test: the same line but with additional (optional) parameters
! and save to an external file
!----------------------------------------------------------------------!
   alpha = 2.0_dp
   mHT = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,  &
      Sw,Ylm,Xlm,alpha)
   write(*,'(A,22X,2F22.15)') "The same function with optional "       &
      // "parameters: ", real(mHT), aimag(mHT)
!----------------------------------------------------------------------!
! Generate the mHT profile from -5 GammaD to +5 GammaD
!----------------------------------------------------------------------!
   open(unit=unit_, file='mHT_profile_H2He.txt', status='replace',     &
      action='write')

   write(unit_,'(A)')'# Frequency       Real(mHT)             Imag(mHT)'

   do point = 1, 1001
      mHT_tabulated(point) = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2, &
         NuOptRe,NuOptIm,nu_tabulated(point),Sw,Ylm,Xlm,alpha)

      write(unit_, '(F15.8, 2F22.15)') nu_tabulated(point),            &
         real(mHT_tabulated(point)), aimag(mHT_tabulated(point))

   end do

   close(unit_)
!----------------------------------------------------------------------!
end program test_mHT
