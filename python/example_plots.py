try: 
  from mHT import mHTprofile
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nexample_plots.py: mHT not found. Make sure catalog with mHT module is in same catalog as running script and is named "mHT"!')
try: 
  import numpy as np
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nexample_plots.py: Numpy not found. Numpy module is needed to run this example.')
try: 
  from matplotlib import pyplot as plt
except ImportError as msg: 
  raise SystemExit (str(msg) + '\nexample_plots.py: Matplotlib not found. Matplotlib module is needed to run this example.')

# Example parameters of the S(1) 3-0 line of H2 perturbed by Ar (Source: 10.1063/5.0139229)
nu0        = 112265.5949 # Unperturbed line position in cm-1.
GammaD     = 35.1e-3     # Doppler broadening in cm-1
Gamma0_Ar  = 11.3e-3     # Speed-averaged line-width in cm-1.
Gamma2_Ar  = 0.374e-3    # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_Ar  = -26.4e-3    # Speed-averaged line-shift in cm-1.
Delta2_Ar  = 17.8e-3     # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_Ar = 72.1e-3     # Real part of the Dicke parameter in cm-1.
NuOptIm_Ar = -16.1e-3    # Imaginary part of the Dicke parameter in cm-1.
nu         = nu0+1       # Current wavenumber in cm-1.

# Optional parameters
Xlm      = 1.0e-3 # Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
Ylm      = 0.5e-3 # Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
alpha_Ar = 20     # Perturber-to-absorber mass ratio, dimensionless.

# Example parameters of the S(1) 3-0 line of H2 perturbed by He (Source: 10.1103/PhysRevA.101.052705)
Gamma0_He  = 11.7e-3  # Speed-averaged line-width in cm-1.
Gamma2_He  = 5.4e-3   # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_He  = 30.5e-3  # Speed-averaged line-shift in cm-1.
Delta2_He  = 12.4e-3  # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_He = 38.0e-3  # Real part of the Dicke parameter in cm-1.
NuOptIm_He = -17.5e-3 # Imaginary part of the Dicke parameter in cm-1.

# Optional parameters
alpha_He = 2 # Perturber-to-absorber mass ratio, dimensionless.

# Preparing tables for plotting 
x        = nu0 + np.linspace(-5*GammaD,+5*GammaD,1001)
y_Ar_Abs = [mHTprofile(nu0,GammaD,Gamma0_Ar,Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu,Ylm,Xlm,alpha_Ar,False) for nu in x]
y_He_Dis = [mHTprofile(nu0,GammaD,Gamma0_He,Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,nu,Ylm,Xlm,alpha_He,True) for nu in x]
y_He_Abs = [mHTprofile(nu0,GammaD,Gamma0_He,Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,nu,Ylm,Xlm,alpha_He,False) for nu in x]

# Plotting absorption and dispersion profile of the H2-He system 
plt.plot(x,y_He_Abs)
plt.plot(x,y_He_Dis)
plt.grid(True)
plt.title("Absorption (blue) and dispersion (orange) part for mHT output for H$_2$-He")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.gca().set_aspect(0.0075)
plt.show()

plt.clf()
# Plotting absorption profiles for H2-Ar and H2-He cases
plt.plot(x,y_Ar_Abs)
plt.plot(x,y_He_Abs)
plt.grid(True)
plt.title("Absorption mHT profiles for H$_2$-Ar (blue) and H$_2$-He (orange).")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.gca().set_aspect(0.0075)
plt.show()
