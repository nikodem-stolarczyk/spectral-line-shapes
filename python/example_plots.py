import mHT 
import numpy as np
from matplotlib import pyplot as plt
#Cosmetics to match the general commenting style.
# Example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
nu0        = 112265.5949; # Unperturbed line position in cm-1.
GammaD     = 35.1e-3;     # Doppler broadening in cm-1
Gamma0_Ar  = 11.3e-3;     # Speed-averaged line-width in cm-1.
Gamma2_Ar  = 0.374e-3;    # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_Ar  = -26.4e-3;    # Speed-averaged line-shift in cm-1.
Delta2_Ar  = 17.8e-3;     # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_Ar = 72.1e-3;     # Real part of the Dicke parameter in cm-1.
NuOptIm_Ar = -16.1e-3;    # Imaginary part of the Dicke parameter in cm-1.
nu         = nu0+1;       # Current wavenumber in cm-1.

# Optional parameters
Ylm      = 1.0e-3; # Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
Xlm      = 0.5e-3; # Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless.
alpha_Ar = 20;     # Perturber-to-absorber mass ratio, dimensionless.

# Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference: 10.1103/PhysRevA.101.052705)
Gamma0_He  = 11.7e-3;  # Speed-averaged line-width in cm-1.
Gamma2_He  = 5.4e-3;   # Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_He  = 30.5e-3;  # Speed-averaged line-shift in cm-1.
Delta2_He  = 12.4e-3;  # Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_He = 38.0e-3;  # Real part of the Dicke parameter in cm-1.
NuOptIm_He = -17.5e-3; # Imaginary part of the Dicke parameter in cm-1.

# Optional parameters
alpha_He = 2; # perturber-to-absorber mass ratio, dimensionless.

# Preparing tables for plotting 
x        = nu0 + np.linspace(-5*GammaD,+5*GammaD,1001)
y_Ar_Abs = [mHT.profile(nu0,GammaD,Gamma0_Ar,Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu,Ylm=Ylm,Xlm=Xlm,alpha=alpha_Ar,disp=False) for nu in x]
y_Ar_Dis = [mHT.profile(nu0,GammaD,Gamma0_Ar,Gamma2_Ar,Delta0_Ar,Delta2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu,Ylm=Ylm,Xlm=Xlm,alpha=alpha_Ar,disp=True) for nu in x]
y_He_Abs = [mHT.profile(nu0,GammaD,Gamma0_He,Gamma2_He,Delta0_He,Delta2_He,NuOptRe_He,NuOptIm_He,nu,Ylm=Ylm,Xlm=Xlm,alpha=alpha_He,disp=False) for nu in x]

# Plotting absorption and dispersion profile of the H2-Ar system 
plt.plot(x,y_Ar_Abs)
plt.plot(x,y_Ar_Dis)
plt.grid(True)
plt.title("Absorption (blue) and dispersion (orange) part for mHT output for H$_2$-Ar")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.gca().set_aspect(0.0075)
plt.show()

plt.clf()
# Comparing absorption profile for H2-Ar and H2-He cases.
plt.plot(x,y_Ar_Abs)
plt.plot(x,y_He_Abs)
plt.grid(True)
plt.title("Absorption mHT profiles for H$_2$-Ar (blue) and H$_2$-He (orange).")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.gca().set_aspect(0.0075)
plt.show()
