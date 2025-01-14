import mHT 

# Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference: 10.1103/PhysRevA.101.052705)
nu0     = 112265.5949 # Unperturbed line position in cm-1.
GammaD  = 35.1e-3     # Doppler broadening in cm-1.
Gamma0  = 11.7e-3     # Speed-averaged line-width in cm-1.
Gamma2  = 5.4e-3      # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0  = 30.5e-3     # Speed-averaged line-shift in cm-1.
Delta2  = 12.4e-3     # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe = 38.0e-3     # Real part of the Dicke parameter in cm-1.
NuOptIm = -17.5e-3    # Imaginary part of the Dicke parameter in cm-1.
nu      = nu0+1       # Current wavenumber in cm-1.

print("The output of the mHT function:")
print("absorbtion:")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu))
print("dispersion:")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,disp=True))

# Optional parameters
newYlm   = 1.0e-3 # Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newXlm   = 0.5e-3 # Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newalpha = 2     # Perturber-to-absorber mass ratio, dimensionless (default: 10.0).

print("The output of the mHT function with optional parameters:")
print("absorbtion:")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,Ylm=newYlm,Xlm=newXlm,alpha=newalpha,disp=False))
print("dispersion:")
print(mHT.profile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,Ylm=newYlm,Xlm=newXlm,alpha=newalpha,disp=True))
# This function accepts 4 optional parameters, provided as dictionary inputs.
# Parameters can be specified in any order, and it is not necessary to provide all of them.
# Optional parameters must be passed last in the function call, using their corresponding keyword (e.g., "Ylm=").
