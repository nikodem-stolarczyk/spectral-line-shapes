      function beta(GammaD, NuOptRe, alpha)
c-----------------------------------------------------------------------
c "beta": Beta-Correction  
c Subroutine to compute beta-correction used for hard-collision
c based line-shape profiles to correct NuOptRe value.
c Applicable up to alpha = 5.0; for higher alpha
c values correction neglected. 
c Source: 10.1016/j.jqsrt.2019.106784
c
c Input/Output Parameters of Routine (Arguments or Common)
c-----------------------------------------------------------------------
c GammaD   : Doppler HWHM in cm-1 (Input) 
c NuOptRe  : Real part of the Dicke parameter in cm-1 (Input).
c alpha    : Mass ratio in the molecule. Applicable up to alpha=5.
c
c The function provides one output:
c-----------------------------------------------------------------------
c (1): Value of the beta correction 
c-----------------------------------------------------------------------
      real*8 :: GammaD, NuOptRe, alpha, beta
c-----------------------------------------------------------------------
c the mass ratio for which the beta correction becomes negligible
c-----------------------------------------------------------------------
      parameter (alpha_max = 5.0d0)
c-----------------------------------------------------------------------
      real*8 :: a, b, c, d
c-----------------------------------------------------------------------
      if (alpha < alpha_max) then
c-----------------------------------------------------------------------
         a =  0.0534d0 + 0.1585d0 * exp(-0.4510d0 * alpha)
         b =  1.9595d0 - 0.1258d0 * alpha + 0.0056d0 * alpha**2.0d0
     &     + 0.0050d0 * alpha**3.0d0
         c = -0.0546d0 + 0.0672d0 * alpha - 0.0125d0 * alpha**2.0d0
     &     + 0.0003d0 * alpha**3.0d0
         d =  0.9466d0 - 0.1585d0 * exp(-0.4510d0 * alpha)
         beta = a * tanh(b * log10(0.5d0 * NuOptRe / GammaD)+c) + d
c-----------------------------------------------------------------------
      else
c-----------------------------------------------------------------------
         beta = 1.0d0
c-----------------------------------------------------------------------
      endif
c-----------------------------------------------------------------------
      end function beta
c-----------------------------------------------------------------------
      function profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,
     &   NuOptIm,nu,Ylm,Xlm,alpha,calculate_dispersion)
c-----------------------------------------------------------------------
c "PROFILE_mHT": modified Hartman Tran profile
c Subroutine to compute the complex normalized spectral shape of an 
c isolated line by the mHT model
c
c Input/Output Parameters of Routine (Arguments or Common)
c-----------------------------------------------------------------------
c nu0       : Unperturbed line position in cm-1 (Input).
c GammaD    : Doppler HWHM in cm-1 (Input)
c Gamma0    : Speed-averaged line-width in cm-1 (Input).       
c Gamma2    : Speed dependence of the line-width in cm-1 (Input).
c Delta0    : Speed-averaged line-shift in cm-1 (Input).
c Delta2    : Speed dependence of the line-shift in cm-1 (Input)   
c NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
c NuOptIm   : Imaginary part of the Dicke parameter in cm-1 (Input).    
c nu        : Current WaveNumber of the Computation in cm-1 (Input).
c Ylm       : Imaginary part of the 1st order (Rosenkranz) line
c             mixing coefficients, dimensionless (Input).
c Xlm       : Real part of the 1st order (Rosenkranz) line mixing
c             coefficients, dimensionless (Input).
c alpha     : Mass ratio in the molecule for calculating
c             beta-correction (Input). Applicable up
c             to alpha=5.
c calculate_dispersion_opt : (Input) false by default:
c             "mHT_profile" returns the real part of the mHT
c             profile (absorption). If true, "mHT_profile" returns
c             the imaginary part of the profile (dispersion).
c
c The function provides one output:
c-----------------------------------------------------------------------
c (1) the normalized spectral shape (cm):
c     - if "calculate_dispersion" is false (default), the function
c       returns the absorption profile.
c     - if "calculate_dispersion" is true, the function returns
c       the dispersion profile.
c-----------------------------------------------------------------------
      include 'constants.inc'
c-----------------------------------------------------------------------
      real*8 :: nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe,
     &   NuOptIm, nu, Ylm, Xlm, alpha, profile
      logical :: calculate_dispersion
      !----------------------------------------------------------------!
      parameter (small_threshold = 3.0d-8)
      !----------------------------------------------------------------!
      complex*16 :: calculated_profile, cpf_accurate
      real*8 :: nuD, nuR
      complex*16 :: c2, c0, LM, X, Y, csqY, z1, z2, w1, w2, wX, A,
     &   X_sqrt, z, w
c-----------------------------------------------------------------------
      nuD = GammaD / sqrt_ln2
      nuR = NuOptRe*beta(GammaD,NuOptRe,alpha)
      c2  = cmplx(Gamma2, Delta2)
      c0  = cmplx(Gamma0, Delta0) - 1.5d0*c2 + nuR
     &    + cmplx(0.0d0, NuOptIm)
      LM  = cmplx(1.0d0 + Xlm, Ylm)
c-----------------------------------------------------------------------
      if ( abs(c2) > numerical_zero ) then
c-----------------------------------------------------------------------
         X    = (cmplx(0.0d0, nu0-nu) + c0) / c2
         Y    = 0.25d0*(nuD/c2)**2.0d0
         csqY = 0.5d0*nuD*cmplx(Gamma2, -Delta2)
     &        /(Gamma2**2.0d0 + Delta2**2.0d0)
         if ( abs( Y )  > abs( X  ) * numerical_zero ) then
c-----------------------------------------------------------------------
            z2 = (X+Y)**0.5d0 + csqY
            if  ( abs(X)  > abs(Y)  * small_threshold ) then
               z1 = z2 - 2.0d0 * csqY
            else
               z1 = (cmplx(0.0d0, nu0-nu) + c0) / nuD    
            endif
            w1 = cpf_accurate(-aimag(z1),real(z1))
            w2 = cpf_accurate(-aimag(z2),real(z2))
            A  = square_root_pi/nuD*(w1-w2)
c-----------------------------------------------------------------------
         else
c-----------------------------------------------------------------------
            X_sqrt = (X)**0.5d0
            if (  abs(X) < numerical_infty ) then
               wX = cpf_accurate(-aimag(X_sqrt),real(X_sqrt))
               A  = 2.0d0*(1.0d0 - square_root_pi*X_sqrt*wX)/c2
            else
               A  = (1.0d0/X - 1.5d0/X**2.0d0)/c2
            endif
c-----------------------------------------------------------------------
         endif
c-----------------------------------------------------------------------
      else
c-----------------------------------------------------------------------
         z = (cmplx(0.0d0, nu0-nu) + c0) / nuD
         w = cpf_accurate(-aimag(z),real(z))
         A = w*square_root_pi/nuD
c-----------------------------------------------------------------------
      endif
c-----------------------------------------------------------------------
      calculated_profile  = LM/pi*A/(1-(nuR + cmplx(0.0d0, NuOptIm))*A)
c-----------------------------------------------------------------------
      if (calculate_dispersion) then
         profile = aimag(calculated_profile)
      else
         profile = real(calculated_profile)
      endif
c-----------------------------------------------------------------------
      end function profile
