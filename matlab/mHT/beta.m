function res = beta(GammaD, NuOptRe, alpha)
    % ---------------------------------------- 
    %    Subroutine to compute beta-correction used for hard-collision based line-shape profiles.
    %    To correct NuOptRe value in the profile. Applicable up to alpha = 5.0, for higher alpha
    %    values correction is neglected.
    %    Reference: 10.1016/j.jqsrt.2019.106784
    %
    %    Standard Input Parameters:
    %    --------------------
    %    GammaD    : Doppler HWHM in cm-1. 
    %    NuOptRe   : Real part of the Dicke parameter in cm-1.
    %    alpha     : Mass ratio in the molecule. Applicable up to alpha=5.0, dimensionless.
    %
    %    The function has one output:
    %    --------------------
    %    (1)       : Value of the beta correction, dimensionless.
    % ----------------------------------------
    if alpha < 5.0 % The mass ratio up to which the beta correction is applicable
        res = (0.0534 + 0.1585*exp(-0.4510*alpha)) *tanh( (1.9595 + alpha*(-0.1258 + alpha*( 0.0056 + alpha*0.0050))) *log10(NuOptRe/GammaD)+ (-0.0546 + alpha*( 0.0672 + alpha*(-0.0125 + alpha*0.0003))) ) + (0.9466 - 0.1585*exp(-0.4510*alpha));
    else
        res = 1.0;
    end
end