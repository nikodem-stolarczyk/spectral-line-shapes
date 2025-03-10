from mHT import mHTprofile

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
print(mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu)) # The mHT function output returns absorption profile by default. To display the dispersive profile see example_dispersion.py.

# Optional parameters
newYlm   = 1.0e-3 # Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newXlm   = 0.5e-3 # Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
newalpha = 2      # Perturber-to-absorber mass ratio, dimensionless (default: 10.0).
newdisp  = True   # Boolean trigger for including dispersion profile in the output (default: False).

print("The output of the mHT function with optional parameters:")
print(mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,newYlm,newXlm,newalpha,newdisp))
# This function accepts 4 optional parameters.
# Parameters are not key-worded, so they need to be specified in specific order.
# To specify Nth parameter all parameters from 0 to N-1 have to be specified as well.
