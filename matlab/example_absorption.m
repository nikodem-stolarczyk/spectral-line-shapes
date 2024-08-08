% Importing global variables used in subfunctions
global e pi rp sln2 num0 numinf
load('mhtconstants.mat')

% Example parameters of the S(1) 3-0 line of H2 perturbed by Ar. 
% Reference: 10.1063/5.0139229.
nu0        = 112265.5949; % Unperturbed line position in cm-1.
GamD       =     35.1e-3; % Doppler broadening in cm-1.
Gam0_Ar    =     11.3e-3; % Speed-averaged line-width in cm-1.
Gam2_Ar    =     37.4e-5; % Speed dependence of the line-width in cm-1.
Shift0_Ar  =    -26.4e-3; % Speed-averaged line-shift in cm-1.
Shift2_Ar  =     17.8e-3; % Speed dependence of the line-shift in cm-1.
NuOptRe_Ar =     72.1e-3; % Real part of the Dicke parameter in cm-1.
NuOptIm_Ar =    -16.1e-3; % Imaginary part of the Dicke parameter in cm-1.

% Example frequency point:
nu = nu0 + 1; % Current wavenumber of the computation in cm-1.

% Displaying mHT function output:
[outRe, outIm] = profile(nu0, GamD, Gam0_Ar, Gam2_Ar, Shift0_Ar, Shift2_Ar, NuOptRe_Ar, NuOptIm_Ar, nu);
disp('The output of the mHT function - absorption:');
disp([outRe, outIm]);
