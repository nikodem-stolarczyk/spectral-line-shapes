import mHT 
from mHT.CPF import cpf_fast, cpf_accurate
from matplotlib import pyplot as plt

x=1; #dimensionless
y=1; #dimensionless
print("the output of the cpf_accurate function")
print(cpf_accurate(x,y))

print("the output of the cpf_fast function")
print(cpf_fast(x,y))

#example parameters of the S(1) 3-0 line od H2 perturbed by Ar (reference 10.1063/5.0139229)
nu0 = 112265.5949;#cm-1
GamD = 35.1e-3;#cm-1
Gam0_Ar = 11.3e-3;#cm-1
Gam2_Ar = 0.374e-3;#cm-1
Shift0_Ar = -26.4e-3;#cm-1
Shift2_Ar = 17.8e-3;#cm-1
NuOptRe_Ar = 72.1e-3;#cm-1
NuOptIm_Ar = -16.1e-3#cm-1
nu=nu0+1;#cm-1

print("the output of the mHT function")
print(mHT.profile(nu0,GamD,Gam0_Ar,Gam2_Ar,Shift0_Ar,Shift2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu))

#optional parameters
Sw = 1.0; #dimensionless
Ylm = 1.0e-3; #dimensionless
Xlm = 0.5e-3; #dimensionless
alpha_Ar = 20; #dimensionless

print("the output of the mHT function with optional parameters")
print(mHT.profile(nu0,GamD,Gam0_Ar,Gam2_Ar,Shift0_Ar,Shift2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu,Sw,Ylm,Xlm,alpha_Ar))

#example parameters of the S(1) 3-0 line od H2 perturbed by He (reference 10.1103/PhysRevA.101.052705)
Gam0_He = 11.7e-3;#cm-1
Gam2_He = 5.4e-3;#cm-1
Shift0_He = 30.5e-3;#cm-1
Shift2_He = 12.4e-3;#cm-1
NuOptRe_He = 38.0e-3;#cm-1
NuOptIm_He = -17.5e-3#cm-1

print("the output of the mHT function")
print(mHT.profile(nu0,GamD,Gam0_He,Gam2_He,Shift0_He,Shift2_He,NuOptRe_He,NuOptIm_He,nu))

#optional parameters
alpha_He = 2; #dimensionless

print("the output of the mHT function with optional parameters")
print(mHT.profile(nu0,GamD,Gam0_He,Gam2_He,Shift0_He,Shift2_He,NuOptRe_He,NuOptIm_He,nu,Sw,Ylm,Xlm,alpha_He))

# Preparing tables for plotting
x    = nu0 + np_linspace(-5*GamD,+5*GamD,1001)
y_Ar = [mHT.profile(nu0,GamD,Gam0_Ar,Gam2_Ar,Shift0_Ar,Shift2_Ar,NuOptRe_Ar,NuOptIm_Ar,nu,Sw,Ylm,Xlm,alpha_Ar) for nu in x]
y_He = [mHT.profile(nu0,GamD,Gam0_He,Gam2_He,Shift0_He,Shift2_He,NuOptRe_He,NuOptIm_He,nu,Sw,Ylm,Xlm,alpha_He) for nu in x]

# Plotting complex output
plt.plot(x,list(map(list,zip(*y_Ar)))[0])
plt.plot(x,list(map(list,zip(*y_Ar)))[1])
plt.grid(True)
plt.title("Real (blue) and imaginary (orange) part for mHT output for H$_2$-Ar")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.show()

plt.clf()
# Comparing real parts for H2-Ar and H2-He cases.
plt.plot(x,list(map(list,zip(*y_Ar)))[0])
plt.plot(x,list(map(list,zip(*y_He)))[0])
plt.grid(True)
plt.title("Real part of mHT output for H$_2$-Ar (blue) and H$_2$-He (orange).")
plt.xlabel("Frequency [cm^(-1)]")
plt.ylabel("mHT")
plt.show()
