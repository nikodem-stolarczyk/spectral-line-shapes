% Add subfolders to PATH variable (can be deleted if package is already in
% folder within PATH or files are in the same catalog as this script).
addpath(genpath(join([fileparts(which(mfilename)),"\mHT Package"],"")));

% Specifing display format output to double precision
format long;

% Example parameters for the CPF functions.
x = 1; % Dimensionless
y = 1; % Dimensionless

% Displaying the cpf function outputs
disp('The output of the cpf_accurate function:');
disp(cpf_accurate(x,y));
disp('The output of the cpf_fast function:');
disp(cpf_fast(x,y));
