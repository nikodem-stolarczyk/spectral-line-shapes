function res = beta(GamD, NuOptRe, alpha)
    % "beta": Beta-Correction  
    % Subroutine to compute beta-correction used for hard-collision based line-shape profiles.
    % To correct NuOptRe value in the profile. Applicable up to alpha = 5.0, for higher alpha
    % values correction is neglected.
    % Source: 10.1016/j.jqsrt.2019.106784
    %
    % Input/Output Parameters of Routine (Arguments or Common)
    % GamD      : Doppler HWHM in cm-1 (Input) 
    % NuOptRe   : Real part of the Dicke parameter in cm-1 (Input)
    % alpha     : Mass ratio in the molecule. Applicable up to alpha=5.
    %
    % The function has one output:
    % (1): Value of the beta correction 

    max_alpha = 5.0; % the mass ratio up to which the beta correction is applicable
    if alpha<max_alpha
        a   = 0.0534+0.1585*exp(-0.4510*alpha);
        b   = 1.9595-0.1258*alpha+0.0056*alpha^2+0.0050*alpha^3;
        c   = -0.0546+0.0672*alpha-0.0125*alpha^2+0.0003*alpha^3;
        d   = 0.9466-0.1585*exp(-0.4510*alpha);
        res = a*tanh(b*log10(0.5*NuOptRe/GamD)+c)+d;
    else
        res = 1.0;
    end
end
