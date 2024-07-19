import mHT 
from mHT.CPF import cpf_fast, cpf_accurate
x=1; #dimensionless
y=1; #dimensionless
print("the output of the cpf_accurate function")
print cpf_accurate(x,y)

print("the output of the cpf_fast function")
print cpf_fast(x,y)

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
print(mHT.profile(nu0,GamD,Gam0_Ar,Gam2_Ar,Shift0_Ar,Shift2_Ar,NuOptR_Are,NuOptIm_Ar,nu))

#optional parameters
Sw = 1.0; #dimensionless
Ylm = 1.0e-3; #dimensionless
Xlm = 0.5e-3; #dimensionless
alpha_Ar = 20; #dimensionless

print("the output of the mHT function with optional parameters")
print(mHT.profile(nu0,GamD,Gam0_Ar,Gam2_Ar,Shift0_Ar,Shift2_Ar,NuOptR_Are,NuOptIm_Ar,nu,Sw,Ylm,Xlm,alpha_Ar))

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
