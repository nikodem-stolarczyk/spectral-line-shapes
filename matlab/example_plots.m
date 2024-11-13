% Add subfolders to PATH variable (can be deleted if package is already in
% folder within PATH or files are in the same catalog as this script).
addpath(fullfile(fileparts(which(mfilename)),"mHT_Package"));

% Example parameters of the S(1) 3-0 line of H2 perturbed by Ar. 
% Reference: 10.1063/5.0139229.
nu0        = 112265.5949; % Unperturbed line position in cm-1.
GamD       =     35.1e-3; % Doppler broadening in cm-1.
Gamma0_Ar  =     11.3e-3; % Speed-averaged line-width in cm-1.
Gamma2_Ar  =     37.4e-5; % Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_Ar  =    -26.4e-3; % Speed-averaged line-shift in cm-1.
Delta2_Ar  =     17.8e-3; % Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_Ar =     72.1e-3; % Real part of the Dicke parameter in cm-1.
NuOptIm_Ar =    -16.1e-3; % Imaginary part of the Dicke parameter in cm-1.

% Example parameters of the S(1) 3-0 line of H2 perturbed by He.
% Reference: 10.1103/PhysRevA.101.052705.
Gamma0_He  =     11.7e-3; % Speed-averaged line-width in cm-1.
Gamma2_He  =      5.4e-3; % Quadratic speed dependence parameter of the line-width in cm-1.
Delta0_He  =     30.5e-3; % Speed-averaged line-shift in cm-1.
Delta2_He  =     12.4e-3; % Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe_He =     38.0e-3; % Real part of the Dicke parameter in cm-1.
NuOptIm_He =    -17.5e-3; % Imaginary part of the Dicke parameter in cm-1.

% Preparing tables for plotting
x  = -5*GamD + nu0:GamD/100:5*GamD + nu0;
yR_Ar = zeros(1,length(x));
yI_Ar = zeros(1,length(x));
yR_He = zeros(1,length(x));
yI_He = zeros(1,length(x));
for i=1:length(x)
    [yR_Ar(i), yI_Ar(i)] = profile(nu0, GamD, Gamma0_Ar, Gamma2_Ar, Delta0_Ar, Delta2_Ar, NuOptRe_Ar, NuOptIm_Ar, x(i));
    [yR_He(i), yI_He(i)] = profile(nu0, GamD, Gamma0_He, Gamma2_He, Delta0_He, Delta2_He, NuOptRe_He, NuOptIm_He, x(i));
end

% Plotting primary example
fig1 = figure(1);
fig1.Position = [100 100 640 340];
plot(x,yR_Ar)
hold on
plot(x,yI_Ar)
hold off
title("Primary example","Absorption (blue) and dispersion (orange) part for mHT output for H_2-Ar")
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
title("Secondary example","Absorption mHT profiles for H_2-Ar (blue) and H_2-He (orange).")
xlabel("Frequency [cm^-^1]")
ylabel("mHT")
grid on
