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
Gam0 = 11.3e-3;#cm-1
Gam2 = 0.374e-3;#cm-1
Shift0 = -26.4e-3;#cm-1
Shift2 = 17.8e-3;#cm-1
NuOptRe = 72.1e-3;#cm-1
NuOptIm = -16.1e-3#cm-1
nu=nu0+1;#cm-1

print("the output of the mHT function")
print(mHT.profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu))

#optional parameters
Sw = 1.0; #dimensionless
Ylm = 1.0e-3; #dimensionless
Xlm = 0.5e-3; #dimensionless
alpha = 20; #dimensionless

print("the output of the mHT function with optional parameters")
print(mHT.profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,nu,Sw,Ylm,Xlm,alpha))
