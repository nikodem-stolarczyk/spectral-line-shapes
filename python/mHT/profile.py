from mHT.CPF import cpf_accurate as cpf
from math import log10, tanh
e      = 2.718281828459045   # Euler constant
pi     = 3.141592653589793   # Pi number
rp     = 1.772453850905516   # root square of pi
sln2   = 0.8325546111576977  # root square of natural logarithm of 2
num0   = 1.0e-15             # numerical zero  
numinf = 4.0e3               # numerical infinity 

def beta(GammaD,NuOptRe,alpha):
    """
    #-------------------------------------------------
    #      "beta": Beta-Correction  
    #      Subroutine to compute beta-correction used for hard-collision based line-shape profiles
    #      To correct NuOptRe value in the profile . Applicable up to alpha = 5.0, for higher alpha
    #      values correction neglected. 
    #      Source: 10.1016/j.jqsrt.2019.106784
    #
    #      Input Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      GammaD    : Doppler broadening in cm-1. 
    #      NuOptRe   : Real part of the Dicke parameter in cm-1.
    #      alpha     : Mass ratio in the molecule, applicable up to alpha=5, dimensionless.
    #
    #      The function has one output:
    #      -----------------
    #      (1): Value of the beta correction 
    #-------------------------------------------------
    """
    max_alpha = 5.0 # the mass ratio up to which the beta correction is applicable
    if alpha < max_alpha:
        a = 0.0534 + 0.1585*e**(-0.4510*alpha)
        b = 1.9595 - 0.1258*alpha + 0.0056*alpha**2 + 0.0050*alpha**3
        c =-0.0546 + 0.0672*alpha - 0.0125*alpha**2 + 0.0003*alpha**3
        d = 0.9466 - 0.1585*e**(-0.4510*alpha)
        return a*tanh(b*log10(0.5*NuOptRe/GammaD)+c)+d
    else:
        return 1.0

def profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,**kwargs):
    """
    #-------------------------------------------------
    #      "PROFILE_mHT": modified Hartman Tran profile
    #      Subroutine to compute the complex normalized spectral-line shape using mHT model
    #
    #      Standard Input Parameters
    #      ---------------------------------
    #      nu0       : Unperturbed line position in cm-1.
    #      GammaD    : Doppler broadening in cm-1.
    #      Gamma0    : Speed-averaged line-width in cm-1.       
    #      Gamma2    : Speed dependence of the line-width in cm-1.
    #      Delta0    : Speed-averaged line-shift in cm-1.
    #      Delta2    : Speed dependence of the line-shift in cm-1.   
    #      NuOptRe   : Real part of the complex Dicke parameter in cm-1.
    #      NuOptIm   : Imaginary part of the complex Dicke parameter in cm-1.    
    #      nu        : Current WaveNumber in cm-1.
    #
    #      Optional Dictionary Input Parameters
    #      ---------------------------------
    #      Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5 (default: 10.0).
    #      disp      : Boolean trigger for including dispersion profile in the output (default: False).
    #
    #      The function has two outputs, from which one is optional:
    #      -----------------
    #      (1):        Real part of the normalized spectral shape (cm)
    #      (Optional): Imaginary part of the normalized spectral shape (cm)
    #
    #-------------------------------------------------
    """
    # Optional input parameters definition
    Xlm   = 0.0   if kwargs.get('Xlm')   == None or not(isinstance(kwargs.get('Xlm'),(float,int)))   else float(kwargs.get('Xlm'))
    Ylm   = 0.0   if kwargs.get('Ylm')   == None or not(isinstance(kwargs.get('Ylm'),(float,int)))   else float(kwargs.get('Ylm'))
    alpha = 10.0  if kwargs.get('alpha') == None or not(isinstance(kwargs.get('alpha'),(float,int))) else float(kwargs.get('alpha'))
    disp  = False if not(kwargs.get('disp',False)==True) else True
    
    nuD = GammaD/sln2
    nuR = NuOptRe*beta(GammaD,NuOptRe,alpha)
    c2  = Gamma2 + Delta2*1j
    c0  = Gamma0 + Delta0*1j - 1.5*c2 + nuR + NuOptIm*1j
    LM  = 1 + Xlm + Ylm*1j
    
    if abs(c2) > num0:
        X    = ((nu0-nu)*1j + c0) / c2
        Y    = 0.25*(nuD/c2)**2
        csqY = 0.50*nuD*(Gamma2 - Delta2*1j)/(Gamma2**2 + Delta2**2)
        if abs(Y)>abs(X)*num0:
            z2 = (X+Y)**0.5 + csqY   
            z1 = z2 - 2*csqY if abs(X)>abs(Y)*3e-8 else ((nu0-nu)*1j + c0) / nuD    
            w1 = cpf(-z1.imag,z1.real)
            w2 = cpf(-z2.imag,z2.real)
            A  = rp/nuD*(w1-w2)
        else:
            if abs(X**0.5) < numinf:
                rX = X**0.5
                wX = cpf(-rX.imag,rX.real)
                A  = 2*(1 - rp*rX*wX)/c2
            else:
                A  = (1/X - 1.5/X**2)/c2
    else:
        z = ((nu0-nu)*1j + c0) / nuD
        w = cpf(-z.imag,z.real)
        A = w*rp/nuD
    I = LM/pi*A/(1-(nuR + NuOptIm*1j)*A)
    return (I.real) if not(disp) else (I.real, I.imag)
