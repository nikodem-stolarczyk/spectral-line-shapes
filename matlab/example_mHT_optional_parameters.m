% Add subfolders to PATH variable (can be deleted if package is already in
% folder within PATH or files are in the same catalog as this script).
addpath(fullfile(fileparts(which(mfilename)),"mHT"));

% Specifing display format output to double precision
format long;

% Example parameters of the S(1) 3-0 line of H2 perturbed by He.
% Reference: 10.1103/PhysRevA.101.052705.
nu0     = 112265.5949; % Unperturbed line position in cm-1.
GammaD  =     35.1e-3; % Doppler broadening in cm-1.
Gamma0  =     11.7e-3; % Speed-averaged line-width in cm-1.
Gamma2  =      5.4e-3; % Quadratic speed dependence parameter of the line-width in cm-1.
Delta0  =     30.5e-3; % Speed-averaged line-shift in cm-1.
Delta2  =     12.4e-3; % Quadratic speed dependence parameter of the line-shift in cm-1.
NuOptRe =     38.0e-3; % Real part of the Dicke parameter in cm-1.
NuOptIm =    -17.5e-3; % Imaginary part of the Dicke parameter in cm-1.

% Example frequency point:
nu = nu0 + 1; % Current wavenumber of the computation in cm-1.

% Optional parameters
Ylm   = 1.0e-3; % Dimensionless 
Xlm   = 0.5e-3; % Dimensionless
alpha = 2 ;     % Dimensionless

% Displaying mHT function output with no optional parameters
[outRe, outIm]  = profile(nu0, GamD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu);
disp('The output of the mHT function with optional parameters:');
disp([outRe, outIm]);

% Displaying mHT function output with Ylm optional parameter only
[outRe, outIm]  = profile(nu0, GamD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, Ylm);
disp('The output of the mHT function with optional parameters (Ylm):');
disp([outRe, outIm]);

% Displaying mHT function output with alpha optional parameter only
[outRe, outIm]  = profile(nu0, GamD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, 0.0, 0.0, alpha);
disp('The output of the mHT function with optional parameters (alpha):');
disp([outRe, outIm]);

% Displaying mHT function output with optional parameters
[outRe, outIm]  = profile(nu0, GamD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, Ylm, Xlm, alpha);
disp('The output of the mHT function with optional parameters (all):');
disp([outRe, outIm]);
