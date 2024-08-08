% Importing global variables used in subfunctions
global e pi rp sln2 num0 numinf
load('mhtconstants.mat')

% Example parameters for the CPF functions.
x = 1; % Dimensionless
y = 1; % Dimensionless

% Specifing format output to double precision (display related)
format long;

% Displaying the cpf function outputs
disp('The output of the cpf_accurate function:');
disp(cpf_accurate(x,y));
disp('The output of the cpf_fast function:');
disp(cpf_fast(x,y));
