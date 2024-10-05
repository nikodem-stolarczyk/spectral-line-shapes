import mHT 

#example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
nu0 = 112265.5949;# Unperturbed line position in cm-1.
GammaD = 35.1e-3;# Doppler broadening in cm-1
Gamma0 = 11.3e-3;# Speed-averaged line-width in cm-1.
Gamma2 = 0.374e-3;# Speed dependence of the line-width in cm-1.
Delta0 = -26.4e-3;# Speed-averaged line-shift in cm-1.
Delta2 = 17.8e-3;# Speed dependence of the line-shift in cm-1.
NuOptRe = 72.1e-3;# Real part of the Dicke parameter in cm-1.
NuOptIm = -16.1e-3# Imaginary part of the Dicke parameter in cm-1.
nu=nu0+1;# Current wavenumber of the computation in cm-1.

print("the output of the mHT function")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu)) # The mHT function output is a table of its real and imaginary part, representing the absorption and dispersion profiles, respectively

#optional parameters
newYlm = 1.0e-3; # Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newXlm = 0.5e-3; # Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newalpha = 20; # perturber-to-absorber mass ratio, dimensionless (default: 10.0).
newdisp = True; # Boolean trigger for including dispersion profile in the output (default: False).

print("the output of the mHT function with optional parameters")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,Ylm=newYlm,Xlm=newXlm,alpha=newalpha,disp=newdisp))
# The function is designed with 4 optional parameters called as a dictionary inputs.
# The parameters doesn't have to be called in same order and doesn't have to be called all at once, only one can be called if that's needed.
# Optional parameters has to be called as last parameters in function and have to be called preceded by its keyworkd (e.g. "Ylm=").
