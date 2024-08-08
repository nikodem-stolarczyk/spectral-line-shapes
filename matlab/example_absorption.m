% Add subfolders to PATH variable (can be deleted if package is already in
% folder within PATH or files are in the same catalog as this script).
addpath(genpath(join([fileparts(which(mfilename)),"\mHT Package"],"")));

% Importing global variables used in subfunctions
global e pi rp sln2 num0 numinf
load('mhtconstants.mat')

% Specifing format output to double precision
format long;

% Example parameters of the S(1) 3-0 line of H2 perturbed by Ar. 
% Reference: 10.1063/5.0139229.
nu0     = 112265.5949; % Unperturbed line position in cm-1.
GamD    =     35.1e-3; % Doppler broadening in cm-1.
Gam0    =     11.3e-3; % Speed-averaged line-width in cm-1.
Gam2    =     37.4e-5; % Speed dependence of the line-width in cm-1.
Shift0  =    -26.4e-3; % Speed-averaged line-shift in cm-1.
Shift2  =     17.8e-3; % Speed dependence of the line-shift in cm-1.
NuOptRe =     72.1e-3; % Real part of the Dicke parameter in cm-1.
NuOptIm =    -16.1e-3; % Imaginary part of the Dicke parameter in cm-1.

% Example frequency point:
nu = nu0 + 1; % Current wavenumber of the computation in cm-1.

% Displaying mHT function output:
outRe = profile(nu0, GamD, Gam0_Ar, Gam2_Ar, Shift0_Ar, Shift2_Ar, NuOptRe_Ar, NuOptIm_Ar, nu);
disp('The output of the mHT function - absorption:');
disp(outRe);

% The mHT function output is a matrix of its real and imaginary part. 
% Read as one parameter output returns only its real part, representing the absorption profile.
