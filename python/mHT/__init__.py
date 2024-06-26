from mHT.CPF import cpf_fast as cpf
from math import log10, tanh
e  = 2.718281828459045
pi = 3.141592653589793
rp = 1.772453850905516

def beta(GamD,NuOptRe,alpha):
    """
    #-------------------------------------------------
    #      "beta": Beta-Correction  
    #      Subroutine to compute beta-correction used for hard-collision based line-shape profiles
    #      To correct NuOptRe value in the profile . Applicable up to alpha = 5.0, for higher alpha
    #      values correction neglected. 
    #      Source: 10.1016/j.jqsrt.2019.106784
    #
    #      Input/Output Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      GamD      : Doppler HWHM in cm-1 (Input) 
    #      NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
    #      alpha     : Mass ratio in the molecule. Applicable up to alpha=5.
    #
    #      The function has one output:
    #      -----------------
    #      (1): Value of the beta correction 
    
    #-------------------------------------------------
    """
    if alpha<5.0:
        a = 0.0534 + 0.1585*e**(-0.4510*alpha)
        b = 1.9595 - 0.1258*alpha + 0.0056*alpha**2 + 0.0050*alpha**3
        c =-0.0546 + 0.0672*alpha - 0.0125*alpha**2 + 0.0003*alpha**3
        d = 0.9466 - 0.1585*e**(-0.4510*alpha)
        return a*tanh(b*log10(0.5*NuOptRe/GamD)+c)+d
    else:
        return 1.0

def profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,Sw=1.0,Ylm=0.0,Xlm=0.0,alpha=10.0):
    """
    #-------------------------------------------------
    #      "PROFILE_mHT": modified Hartman Tran profile
    #      Subroutine to compute the complex normalized spectral shape of an 
    #      isolated line by the mHT model
    #
    #      Input/Output Parameters of Routine (Arguments or Common)
    #      ---------------------------------
    #      nu0       : Unperturbed line position in cm-1 (Input).
    #      GamD      : Doppler HWHM in cm-1 (Input)
    #      Gam0      : Speed-averaged line-width in cm-1 (Input).       
    #      Gam2      : Speed dependence of the line-width in cm-1 (Input).
    #      Shift0    : Speed-averaged line-shift in cm-1 (Input).
    #      Shift2    : Speed dependence of the line-shift in cm-1 (Input)   
    #      NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
    #      NuOptIm   : Imaginary part of the Dicke parameter in cm-1 (Input).    
    #      nu        : Current WaveNumber of the Computation in cm-1 (Input).
    #	   Sw		 : Statistical weight 
    #      Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients in cm-1 (Input)
    #      Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients in cm-1 (Input)
    #      alpha     : Mass ratio in the molecule for calculating beta-correction. Applicable up to alpha=5.
    #
    #      The function has two outputs:
    #      -----------------
    #      (1): Real part of the normalized spectral shape (cm)
    #      (2): Imaginary part of the normalized spectral shape (cm)
    #
    #-------------------------------------------------
    """
    nuD = GamD/0.8325546111576977
    nuR = NuOptRe*beta(GamD,NuOptRe,alpha)
    c2  = Gam2 + Shift2*1j
    c0  = Gam0 + Shift0*1j - 1.5*c2 + nuR + NuOptIm*1j
    LM  = 1 + Xlm + Ylm*1j
    
    if abs(c2) != 0.0:
        X    = ((nu0-nu)*1j + c0) / c2
        Y    = 0.25*(nuD/c2)**2
        csqY = 0.50*nuD*(Gam2 - Shift2*1j)/(Gam2**2 + Shift2**2)
        if abs(Y)>abs(X)*1.0e-15:
            z2 = (X+Y)**0.5 + csqY   
            z1 = z2 - 2*csqY if abs(X)>abs(Y)*3e-8 else ((nu0-nu)*1j + c0) / nuD    
            w1 = cpf(-z1.imag,z1.real)
            w2 = cpf(-z2.imag,z2.real)
            A  = rp/nuD*(w1-w2)
        else:
            if abs(X**0.5) < 4.0e3:
                rX = X**0.5
                wX = cpf(-rX.imag,rX.real)
                A  = 2*(1 - rp*rX*wX)/c2
            else:
                A  = (1/X - 1.5/X**2)/c2
    else:
        z = ((nu0-nu)*1j + c0) / nuD
        w = cpf(-z.imag,z.real)
        A = w*rp/nuD
    I = Sw*LM/pi*A/(1-(nuR + NuOptIm*1j)*A)
    return I.real, I.imag
