from mHT.CPF import cpf_accurate as cpf
from math import log10, tanh
import numpy as np
e      = 2.718281828459045   # Euler constant
pi     = 3.141592653589793   # Pi number
rp     = 1.772453850905516   # root square of pi
sln2   = 0.8325546111576977  # root square of natural logarithm of 2
num0   = 1.0e-15             # numerical zero  
numinf = 4.0e3               # numerical infinity 

def beta(GammaD,NuOptRe,alpha):
    """
    # ----------------------------------------
    #      "BETA": Beta-Correction Function
    #      Subroutine to compute beta-correction used for hard-collision based line-shape profiles
    #      To correct NuOptRe value in the profile . Applicable up to alpha = 5.0, for higher alpha
    #      values correction neglected. 
    #      Source: 10.1016/j.jqsrt.2019.106784
    #
    #      Standard Input Parameters:
    #      --------------------
    #      GammaD    : Doppler broadening in cm-1. 
    #      NuOptRe   : Real part of the Dicke parameter in cm-1.
    #      alpha     : Mass ratio in the molecule, applicable up to alpha=5, dimensionless.
    #
    #      The function has one output:
    #      --------------------
    #      (1)       : Value of the beta correction, dimensionless. 
    # ----------------------------------------
    """
    max_alpha = 5.0 # the mass ratio up to which the beta correction is applicable
    if alpha < max_alpha:
        a = 0.0534 + 0.1585*pow(e,-0.4510*alpha)
        b = 1.9595 + alpha*(-0.1258 + alpha*( 0.0056 + alpha*0.0050))
        c =-0.0546 + alpha*( 0.0672 + alpha*(-0.0125 + alpha*0.0003))
        d = 0.9466 - 0.1585*pow(e,-0.4510*alpha)
        return a*tanh(b*log10(NuOptRe/GammaD)+c)+d
    else:
        return 1.0

def profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,**kwargs):
    """
    # ----------------------------------------
    #      "PROFILE": modified Hartmann Tran profile
    #      Subroutine to compute the complex normalized spectral-line shape using mHT model
    #
    #      Standard Input Parameters:
    #      --------------------
    #      nu0       : Unperturbed line position in cm-1.
    #      GammaD    : Doppler broadening in cm-1.
    #      Gamma0    : Speed-averaged line-width in cm-1.       
    #      Gamma2    : Quadratic speed dependence parameter of the line-width in cm-1.
    #      Delta0    : Speed-averaged line-shift in cm-1.
    #      Delta2    : Quadratic speed dependence parameter of the line-shift in cm-1.   
    #      NuOptRe   : Real part of the complex Dicke parameter in cm-1.
    #      NuOptIm   : Imaginary part of the complex Dicke parameter in cm-1.    
    #      nu        : Current WaveNumber in cm-1.
    #
    #      Optional Input Parameters:
    #      --------------------
    #      Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless (default: 10.0).
    #      disp      : Boolean trigger for including dispersion profile in the output (default: False).
    #
    #      The function has one outputs:
    #      --------------------
    #      (1)       : Real or imaginary (depending on disp value) part of the normalized spectral shape in cm.
    #
    # ----------------------------------------
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
        Y    = 0.25*pow((nuD/c2),2)
        csqY = 0.50*nuD*(Gamma2 - Delta2*1j)/(pow(Gamma2,2) + pow(Delta2,2))
        if abs(Y)>abs(X)*num0:
            z2 = pow((X+Y),0.5) + csqY   
            z1 = z2 - 2*csqY if abs(X)>abs(Y)*3e-8 else ((nu0-nu)*1j + c0) / nuD    
            w1 = cpf(-z1.imag,z1.real)
            w2 = cpf(-z2.imag,z2.real)
            A  = rp/nuD*(w1-w2)
        else:
            if abs(X**0.5) < numinf:
                rX = pow(X,0.5)
                wX = cpf(-rX.imag,rX.real)
                A  = 2*(1 - rp*rX*wX)/c2
            else:
                A  = (1/X - 1.5/X**2)/c2
    else:
        z = ((nu0-nu)*1j + c0) / nuD
        w = cpf(-z.imag,z.real)
        A = w*rp/nuD
    I = LM/pi*A/(1-(nuR + NuOptIm*1j)*A)
    return (I.real) if not(disp) else (I.imag)


def profile_vector(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,**kwargs):
    """
    # ----------------------------------------
    #      "PROFILE": modified Hartmann Tran profile
    #      Subroutine to compute the complex normalized spectral-line shape using mHT model
    #
    #      Standard Input Parameters:
    #      --------------------
    #      nu0       : Unperturbed line position in cm-1.
    #      GammaD    : Doppler broadening in cm-1.
    #      Gamma0    : Speed-averaged line-width in cm-1.       
    #      Gamma2    : Quadratic speed dependence parameter of the line-width in cm-1.
    #      Delta0    : Speed-averaged line-shift in cm-1.
    #      Delta2    : Quadratic speed dependence parameter of the line-shift in cm-1.   
    #      NuOptRe   : Real part of the complex Dicke parameter in cm-1.
    #      NuOptIm   : Imaginary part of the complex Dicke parameter in cm-1.    
    #      nu        : Current WaveNumber in cm-1.
    #
    #      Optional Input Parameters:
    #      --------------------
    #      Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    #      alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless (default: 10.0).
    #      disp      : Boolean trigger for including dispersion profile in the output (default: False).
    #
    #      The function has one outputs:
    #      --------------------
    #      (1)       : Real or imaginary (depending on disp value) part of the normalized spectral shape in cm.
    #
    # ----------------------------------------
    """
    # Optional input parameters definition
    Xlm   = 0.0   if kwargs.get('Xlm')   == None or not(isinstance(kwargs.get('Xlm'),(float,int)))   else float(kwargs.get('Xlm'))
    Ylm   = 0.0   if kwargs.get('Ylm')   == None or not(isinstance(kwargs.get('Ylm'),(float,int)))   else float(kwargs.get('Ylm'))
    alpha = 10.0  if kwargs.get('alpha') == None or not(isinstance(kwargs.get('alpha'),(float,int))) else float(kwargs.get('alpha'))
    disp  = False if not(kwargs.get('disp',False)==True) else True
    
    nu = np.array(nu)
    
    nuD = GammaD/sln2
    nuR = NuOptRe*beta(GammaD,NuOptRe,alpha)
    c2  = Gamma2 + Delta2*1j
    c0  = Gamma0 + Delta0*1j - 1.5*c2 + nuR + NuOptIm*1j
    LM  = 1 + Xlm + Ylm*1j
    
    if abs(c2) > num0:
        X    = ((nu0-nu)*1j + c0) / c2
        Y    = 0.25*pow((nuD/c2),2)
        csqY = 0.50*nuD*(Gamma2 - Delta2*1j)/(pow(Gamma2,2) + pow(Delta2,2))

        condition = np.abs(Y) > np.abs(X) * num0
        z1 = np.where(condition, np.sqrt(X + Y) + csqY, ((nu0 - nu) * 1j + c0) / nuD)
        z2 = np.where(condition, z1 - 2 * csqY, np.sqrt(X + Y) + csqY)
        w1 = cpf(-z1.imag, z1.real)
        w2 = cpf(-z2.imag, z2.real)
        A = rp / nuD * (w1 - w2)

        abs_X_sqrt = np.abs(np.sqrt(X))
        A = np.where(abs_X_sqrt < numinf, 2 * (1 - rp * np.sqrt(X) * cpf(-np.sqrt(X).imag, np.sqrt(X).real)) / c2,(1 / X - 1.5 / X ** 2) / c2)

    else:
        z = ((nu0-nu)*1j + c0) / nuD
        w = cpf(-z.imag,z.real)
        A = w*rp/nuD
    I = LM/pi*A/(1-(nuR + NuOptIm*1j)*A)
    return (I.real) if not(disp) else (I.imag)
