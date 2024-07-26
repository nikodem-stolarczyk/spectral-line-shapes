module spectral_module
!----------------------------------------------------------------------!
! This module provides 2 functions 
! 1. beta - computes the beta-correction used for the hard-collision
!    based line-shape profiles
! 2. profile - computes the modified Hartmann-Tran (mHT) profile
!----------------------------------------------------------------------!
   use, intrinsic :: iso_fortran_env, only: int32, dp => real64
   use cpf_module, only: cpf_accurate
   implicit none
!----------------------------------------------------------------------!
   real(dp), parameter :: e  = 2.718281828459045_dp
   real(dp), parameter :: pi = 3.141592653589793_dp
   real(dp), parameter :: square_root_pi = 1.772453850905516_dp
   real(dp), parameter :: sqrt_ln2 = 0.8325546111576977_dp
   real(dp), parameter :: numerical_zero = 1e-15_dp
   real(dp), parameter :: numerical_infty = 4e3_dp
!----------------------------------------------------------------------!
   contains
!----------------------------------------------------------------------!
   function beta(GamD, NuOptRe, alpha) result(beta_result)
      !----------------------------------------------------------------!
      ! "beta": Beta-Correction  
      ! Subroutine to compute beta-correction used for hard-collision
      ! based line-shape profiles to correct NuOptRe value.
      ! Applicable up to alpha = 5.0; for higher alpha
      ! values correction neglected. 
      ! Source: 10.1016/j.jqsrt.2019.106784
      !
      ! Input/Output Parameters of Routine (Arguments or Common)
      ! ---------------------------------------------------------------!
      ! GamD      : Doppler HWHM in cm-1 (Input) 
      ! NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
      ! alpha     : Mass ratio in the molecule. Applicable up to alpha=5.
      !
      ! The function provides one output:
      ! ---------------------------------------------------------------!
      ! (1): Value of the beta correction 
      !----------------------------------------------------------------!
      real(dp), intent(in) :: GamD, NuOptRe, alpha
      real(dp) :: beta_result
      !----------------------------------------------------------------!
      ! the mass ratio for which the beta correction becomes negligible
      !----------------------------------------------------------------!
      real(dp), parameter :: max_alpha = 5.0_dp
      !----------------------------------------------------------------!   
      real(dp) :: a, b, c, d
      !----------------------------------------------------------------!
      if (alpha < max_alpha) then
         !-------------------------------------------------------------!
         a =  0.0534_dp + 0.1585_dp * exp(-0.4510_dp * alpha)
         b =  1.9595_dp - 0.1258_dp * alpha + 0.0056_dp * alpha**2.0_dp&
           + 0.0050_dp * alpha**3.0_dp
         c = -0.0546_dp + 0.0672_dp * alpha - 0.0125_dp * alpha**2.0_dp&
           + 0.0003_dp * alpha**3.0_dp
         d =  0.9466_dp - 0.1585_dp * exp(-0.4510_dp * alpha)
         beta_result = a * tanh(b * log10(0.5_dp * NuOptRe / GamD) + c)&
           + d
         !-------------------------------------------------------------!
      else
         !-------------------------------------------------------------!
         beta_result = 1.0_dp
         !-------------------------------------------------------------!
      endif
      !----------------------------------------------------------------!
   end function beta
!----------------------------------------------------------------------!
   function profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,&
      Ylm_opt,Xlm_opt,alpha_opt) result(mHT_profile)
      !----------------------------------------------------------------!
      ! "PROFILE_mHT": modified Hartman Tran profile
      ! Subroutine to compute the complex normalized spectral shape of an 
      ! isolated line by the mHT model
      !
      ! Input/Output Parameters of Routine (Arguments or Common)
      ! ---------------------------------------------------------------!
      ! nu0       : Unperturbed line position in cm-1 (Input).
      ! GamD      : Doppler HWHM in cm-1 (Input)
      ! Gam0      : Speed-averaged line-width in cm-1 (Input).       
      ! Gam2      : Speed dependence of the line-width in cm-1 (Input).
      ! Shift0    : Speed-averaged line-shift in cm-1 (Input).
      ! Shift2    : Speed dependence of the line-shift in cm-1 (Input)   
      ! NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
      ! NuOptIm   : Imaginary part of the Dicke parameter in cm-1 (Input).    
      ! nu        : Current WaveNumber of the Computation in cm-1 (Input).
      ! Ylm       : Imaginary part of the 1st order (Rosenkranz) line
      !             mixing coefficients, dimensionless (Input)
      ! Xlm       : Real part of the 1st order (Rosenkranz) line mixing
      !             coefficients, dimensionless (Input)
      ! alpha     : Mass ratio in the molecule for calculating
      !             beta-correction. Applicable up to alpha=5.
      !
      ! The function provides one, complex output:
      ! ---------------------------------------------------------------!
      ! (1) the normalized spectral shape (cm)
      !----------------------------------------------------------------!
      real(dp), intent(in) :: nu0, GamD, Gam0, Gam2, Shift0, Shift2,   &
         NuOptRe, NuOptIm, nu
      real(dp), intent(in), optional :: Ylm_opt, Xlm_opt, alpha_opt
      complex(dp) :: mHT_profile
      !----------------------------------------------------------------!
      real(dp), parameter :: small_threshold = 3e-8_dp
      !----------------------------------------------------------------!
      real(dp) :: Ylm, Xlm, alpha
      real(dp) :: nuD, nuR
      complex(dp) :: c2, c0, LM, X, Y, csqY, z1, z2, w1, w2, wX, A,    &
         X_sqrt, z, w
      !----------------------------------------------------------------!
      if (present(Ylm_opt)) then
         Ylm = Ylm_opt
      else
         Ylm = 0.0_dp
      endif
      !----------------------------------------------------------------!
      if (present(Xlm_opt)) then
         Xlm = Xlm_opt
      else
         Xlm = 0.0_dp
      endif 
      !----------------------------------------------------------------!
      if (present(alpha_opt)) then
         alpha = alpha_opt
      else
         alpha = 10.0_dp
      endif 
      !----------------------------------------------------------------!
      nuD = GamD / sqrt_ln2
      nuR = NuOptRe*beta(GamD,NuOptRe,alpha)
      c2  = cmplx(Gam2, Shift2, kind=dp)
      c0  = cmplx(Gam0, Shift0, kind=dp) - 1.5_dp*c2 + nuR             &
          + cmplx(0.0_dp, NuOptIm, kind=dp)
      LM  = cmplx(1.0_dp + Xlm, Ylm, kind=dp)
      !----------------------------------------------------------------!
      if ( abs(c2) > numerical_zero ) then
         !-------------------------------------------------------------!
         X    = (cmplx(0_dp, nu0-nu, kind=dp) + c0) / c2
         Y    = 0.25_dp*(nuD/c2)**2.0_dp
         csqY = 0.5_dp*nuD*cmplx(Gam2, -Shift2, kind=dp)               &
              /(Gam2**2.0_dp + Shift2**2.0_dp)
         if ( abs( Y )  > abs( X  ) * numerical_zero ) then
            !----------------------------------------------------------!
            z2 = (X+Y)**0.5_dp + csqY
            if  ( abs(X)  > abs(Y)  * small_threshold ) then
               z1 = z2 - 2.0_dp * csqY
            else
               z1 = (cmplx(0.0_dp, nu0-nu, kind=dp) + c0) / nuD    
            endif
            w1 = cpf_accurate(-aimag(z1),real(z1))
            w2 = cpf_accurate(-aimag(z2),real(z2))
            A  = square_root_pi/nuD*(w1-w2)
            !----------------------------------------------------------!
         else
            !----------------------------------------------------------!
            X_sqrt = (X)**0.5_dp
            if (  abs(X) < numerical_infty ) then
               wX = cpf_accurate(-aimag(X_sqrt),real(X_sqrt))
               A  = 2.0_dp*(1.0_dp - square_root_pi*X_sqrt*wX)/c2
            else
               A  = (1.0_dp/X - 1.5_dp/X**2.0_dp)/c2
            endif
            !----------------------------------------------------------!
         endif
         !-------------------------------------------------------------!
      else
         !-------------------------------------------------------------!
         z = (cmplx(0.0_dp, nu0-nu, kind=dp) + c0) / nuD
         w = cpf_accurate(-aimag(z),real(z))
         A = w*square_root_pi/nuD
         !-------------------------------------------------------------!
      endif
      
      mHT_profile = LM/pi*A/(1-(nuR + cmplx(0.0_dp, NuOptIm, kind=dp))*A)

   end function profile 

end module spectral_module
