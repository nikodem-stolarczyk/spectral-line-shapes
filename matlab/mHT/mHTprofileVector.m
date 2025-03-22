function [real_part, imag_part] = mHTprofileVector(nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, varargin)
    % ---------------------------------------- 
    %    Subroutine to compute the complex normalized spectral shape of an 
    %    isolated line by the modified Hartmann Tran (mHT) model
    % 
    %    To decrease code execution time, following 
    %    numeric values were introduced explicitly:
    %     * 1/Pi           = 0.3183098861837907
    %     * Sqrt[Pi]       = 1.772453850905516
    %     * 2 Sqrt[Pi]     = 3.5449077018110318
    %     * 1/Sqrt[Log[2]] = 1.2011224087864498
    %
    %    Standard Input Parameters:
    %    --------------------
    %    nu0       : Unperturbed line position in cm-1.
    %    GammaD    : Doppler HWHM in cm-1.
    %    Gamma0    : Speed-averaged line-width in cm-1.       
    %    Gamma2    : Quadratic speed dependence parameter of the line-width in cm-1.
    %    Delta0    : Speed-averaged line-shift in cm-1.
    %    Delta2    : Quadratic speed dependence parameter of the line-shift in cm-1.   
    %    NuOptRe   : Real part of the Dicke parameter in cm-1.
    %    NuOptIm   : Imaginary part of the Dicke parameter in cm-1.    
    %    nu        : Current wavenumber in cm-1.
    %
    %    Optional Input Parameters:
    %    --------------------
    %    Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    %    Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    %    alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha = 5.0, dimensionless (default: 10.0).
    %
    %    The function has two outputs:
    %    --------------------
    %    real_part : Real part of the normalized spectral shape in cm.
    %    imag_part : Imaginary part of the normalized spectral shape in cm.
    % ---------------------------------------- 
    nuD = 1.2011224087864498*GammaD;

    switch nargin-9 % Number of optional parameters
        case 0 
            LM  = 1.0;
            nuR = NuOptRe;
        case 1
            LM  = 1-varargin{1}*1i;
            nuR = NuOptRe;
        case 2
            LM  = 1+varargin{2}-varargin{1}*1i;
            nuR = NuOptRe;
        otherwise
            LM  = 1+varargin{2}-varargin{1}*1i;
            nuR = NuOptRe*beta(GammaD, NuOptRe, varargin{3});
    end

    c2  = Gamma2+Delta2*1i;
    c0  = Gamma0+Delta0*1i-1.5*c2+nuR+NuOptIm*1i;
    if abs(c2) > 1.0e-9 % Limit where speed dependence impact is lower than numerical noise level
        X    = ((nu0-nu)*1i+c0)/c2;
        Y    = 0.25*(nuD/c2)^2;
        csqY = 0.50*nuD*(Gamma2-Delta2*1i)/(Gamma2^2+Delta2^2);
        if real(c0)>1e-9
            vecabsX = min(abs(X));
        else
            vecabsX = median(abs(X));
        end
        if abs(Y)>vecabsX*1.0e-15 % Numerical zero
            z2 = sqrt(X+Y)+csqY;   
            if vecabsX>abs(Y)*3e-8 %
                z1 = z2-2*csqY;
            else
                z1 = ((nu0-nu)*1i+c0)/nuD;
            end
            w1 = cpf_accurate_vector(-imag(z1), real(z1));
            w2 = cpf_accurate_vector(-imag(z2), real(z2));
            A  = 1.772453850905516/nuD.*(w1-w2);
        else
            if real(c0)>1e-9
                vecabssqX = min(abs(sqrt(X)));
            else
                vecabssqX = median(abs(sqrt(X)));
            end
            if vecabssqX < 4000 % Numerical infinity
                rX = sqrt(X);
                wX = cpf_accurate_vector(-imag(rX), real(rX));
                A  = (2-3.5449077018110318.*rX.*wX)/c2;
            else
                A  = (1-1.5./X)./X/c2;
            end
        end
    else
        z = ((nu0-nu)*1i+c0)/nuD;
        w = cpf_accurate_vector(-imag(z), real(z));
        A = 1.772453850905516.*w/nuD;
    end
    I = 0.3183098861837907*LM.*A./(1-(nuR+NuOptIm*1i).*A);
    
    % Output
    real_part = real(I);
    imag_part = -imag(I);
end
