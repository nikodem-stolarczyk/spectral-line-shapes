import mHT 

#example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference 10.1063/5.0139229)
nu0 = 112265.5949;# Unperturbed line position in cm-1.
GamD = 35.1e-3;# Doppler broadening in cm-1
Gam0 = 11.3e-3;# Speed-averaged line-width in cm-1.
Gam2 = 0.374e-3;# Speed dependence of the line-width in cm-1.
Shift0 = -26.4e-3;# Speed-averaged line-shift in cm-1.
Shift2 = 17.8e-3;# Speed dependence of the line-shift in cm-1.
NuOptRe = 72.1e-3;# Real part of the Dicke parameter in cm-1.
NuOptIm = -16.1e-3# Imaginary part of the Dicke parameter in cm-1.
nu=nu0+1;# Current wavenumber of the computation in cm-1.

print("the output of the mHT function - absorption")
print(mHT.profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu)[0]) # The mHT function output is a table of its real and imaginary part. The part [0] is its real part, representing the absorption profile.
