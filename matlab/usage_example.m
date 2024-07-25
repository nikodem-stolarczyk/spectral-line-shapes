% Importing global variables used in subfunctions
global e pi rp sln2 num0 numinf
load('mhtconstants.mat')

% Example parameters for the CPF functions.
x = 1; % Dimensionless
y = 1; % Dimensionless

% Specifing format output to double precision
format long;

% Displaying the cpf function outputs
disp('The output of the cpf_accurate function:');
disp(cpf_accurate(x,y));
disp('The output of the cpf_fast function:');
disp(cpf_fast(x,y));

% Example parameters of the S(1) 3-0 line of H2 perturbed by Ar (reference: 10.1063/5.0139229) 
nu0        = 112265.5949; % Unperturbed line position in cm-1
GamD       =     35.1e-3; % Doppler broadening in cm-1
Gam0_Ar    =     11.3e-3; % Speed-averaged line-width in cm-1
Gam2_Ar    =     37.4e-5; % Speed dependence of the line-width in cm-1
Shift0_Ar  =    -26.4e-3; % Speed-averaged line-shift in cm-1
Shift2_Ar  =     17.8e-3; % Speed dependence of the line-shift in cm-1
NuOptRe_Ar =     72.1e-3; % Real part of the Dicke parameter in cm-1
NuOptIm_Ar =    -16.1e-3; % Imaginary part of the Dicke parameter in cm-1

% Example frequency point:
nu = nu0 + 1; % Current wavenumber of the computation in cm-1

% Displaying mHT function output
[outRe, outIm] = profile(nu0, GamD, Gam0_Ar, Gam2_Ar, Shift0_Ar, Shift2_Ar, NuOptRe_Ar, NuOptIm_Ar, nu);
disp('The output of the mHT function:');
disp([outRe, outIm]);

% Optional parameters
Ylm      = 1.0e-3; % Dimensionless 
Xlm      = 0.5e-3; % Dimensionless
alpha_Ar = 20;     % Dimensionless

% Displaying mHT function output with opt parameters
[outRe, outIm] = profile(nu0, GamD, Gam0_Ar, Gam2_Ar, Shift0_Ar, Shift2_Ar, NuOptRe_Ar, NuOptIm_Ar, nu, Ylm, Xlm, alpha_Ar);
disp('The output of the mHT function with optional parameters:');
disp([outRe, outIm]);

% Example parameters of the S(1) 3-0 line of H2 perturbed by He (reference 10.1103/PhysRevA.101.052705)
Gam0_He    =  11.7e-3; % cm-1
Gam2_He    =   5.4e-3; % cm-1
Shift0_He  =  30.5e-3; % cm-1
Shift2_He  =  12.4e-3; % cm-1
NuOptRe_He =  38.0e-3; % cm-1
NuOptIm_He = -17.5e-3; % cm-1

% Displaying mHT function output
[outRe, outIm] = profile(nu0, GamD, Gam0_He, Gam2_He, Shift0_He, Shift2_He, NuOptRe_He, NuOptIm_He, nu);
disp('The real-part output of the mHT function (second example):');
disp([outRe, outIm]);

% Optional parameters
alpha_He = 2; % Dimensionless

% Displaying mHT function output with opt parameters
[outRe, outIm] = profile(nu0, GamD, Gam0_He, Gam2_He, Shift0_He, Shift2_He, NuOptRe_He, NuOptIm_He, nu, Ylm, Xlm, alpha_He);
disp('The output of the mHT function with optional parameters (second example):');
disp([outRe, outIm]);

% Preparing tables for plotting
x  = -5*GamD + nu0:GamD/100:5*GamD + nu0;
yR_Ar = zeros(1,length(x));
yI_Ar = zeros(1,length(x));
yR_He = zeros(1,length(x));
yI_He = zeros(1,length(x));
for i=1:length(x)
    [yR_Ar(i), yI_Ar(i)] = profile(nu0, GamD, Gam0_Ar, Gam2_Ar, Shift0_Ar, Shift2_Ar, NuOptRe_Ar, NuOptIm_Ar, x(i));
    [yR_He(i), yI_He(i)] = profile(nu0, GamD, Gam0_He, Gam2_He, Shift0_He, Shift2_He, NuOptRe_He, NuOptIm_He, x(i));
end

% Plotting primary example
fig1 = figure(1);
fig1.Position = [100 100 640 340];
plot(x,yR_Ar)
hold on
plot(x,yI_Ar)
hold off
title("Primary example","Real (blue) and imaginary (orange) part of the mHT output")
xlabel("Frequency [cm^-^1]")
ylabel("mHT")
grid on

% Plotting secondary example
fig2 = figure(2);
fig2.Position = [750 100 640 340];
plot(x,yR_Ar)
hold on
plot(x,yR_He)
hold off
title("Secondary example","Real part of the H_2-Ar (blue) and the H_2-He (orange) mHT output")
xlabel("Frequency [cm^-^1]")
ylabel("mHT")
grid on
