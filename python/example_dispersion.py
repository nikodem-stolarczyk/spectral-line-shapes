try: 
  from mHT import mHTprofile
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nexample_dispersion.py: mHT not found. Make sure catalog with mHT module is in same catalog as running script and is named "mHT"!')

# Example parameters of the S(1) 3-0 line of H2 perturbed by He (Source: 10.1103/PhysRevA.101.052705)
nu0     = 112265.5949 # Unperturbed line position in cm-1.
GammaD  = 35.1e-3     # Doppler broadening in cm-1.
Gamma0  = 11.7e-3     # Speed-averaged line-width in cm-1.
Gamma2  = 5.4e-3      # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0  = 30.5e-3     # Speed-averaged line-shift in cm-1.
Delta2  = 12.4e-3     # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe = 38.0e-3     # Real part of the Dicke parameter in cm-1.
NuOptIm = -17.5e-3    # Imaginary part of the Dicke parameter in cm-1.
nu      = nu0+1       # Current wavenumber in cm-1.

print("The output of the mHT function - dispersion:")
print(mHTprofile(nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,0,0,10,True))
