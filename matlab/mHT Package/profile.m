function [real_part, imag_part] = profile(nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, Ylm, Xlm, alpha)
    % ---------------------------------------- 
    %    Subroutine to compute the complex normalized spectral shape of an 
    %    isolated line by the modified Hartman Tran (mHT) model
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
    %    alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless (default: 10.0).
    %
    %    The function has two outputs:
    %    --------------------
    %    real_part : Real part of the normalized spectral shape in cm.
    %    imag_part : Imaginary part of the normalized spectral shape in cm.
    % ---------------------------------------- 
    arguments (Input)
        nu0     (1,1) double
        GammaD  (1,1) double
        Gamma0  (1,1) double
        Gamma2  (1,1) double
        Delta0  (1,1) double
        Delta2  (1,1) double
        NuOptRe (1,1) double
        NuOptIm (1,1) double
        nu      (1,1) double
        Ylm     (1,1) double = 0.0
        Xlm     (1,1) double = 0.0
        alpha   (1,1) double = 10.0
    end    
    arguments (Output)
        real_part (1,1) double
        imag_part (1,1) double
    end
    
    % Define mathematical, numerical and function handles variables
    pi      = 3.141592653589793;  % Pi number
    rp      = 1.772453850905516;  % Root square of pi
    sln2    = 0.8325546111576977; % Root square of natural logarithm of 2
    num0    = 1.0e-15;            % Numerical zero
    numinf  = 4000;               % Numerical infinity
    cpf     = @cpf_accurate;      % CPF choice

    % Ternary operator
    ternary = @(varargin) varargin{length(varargin)-varargin{1}}; 
    
    % The function
    nuD = GammaD / sln2 ;
    nuR = NuOptRe*beta(GammaD, NuOptRe, alpha);
    c2  = Gamma2+Delta2*1i;
    c0  = Gamma0+Delta0*1i-1.5*c2+nuR+NuOptIm*1i;
    LM  = 1+Xlm+Ylm*1i;
    if abs(c2) ~= num0
        X    = ((nu0-nu)*1i+c0)/c2;
        Y    = 0.25*(nuD/c2)^2;
        csqY = 0.50*nuD*(Gamma2-Delta2*1i)/(Gamma2^2+Delta2^2);
        if abs(Y)>abs(X)*num0
            z2 = sqrt(X+Y)+csqY;   
            z1 = ternary(abs(X)>abs(Y)*3e-8,z2-2*csqY, ((nu0-nu)*1i+c0)/nuD);
            w1 = cpf(-imag(z1), real(z1));
            w2 = cpf(-imag(z2), real(z2));
            A  = rp/nuD*(w1-w2);
        else
            if abs(sqrt(X)) < numinf
                rX = sqrt(X);
                wX = cpf(-imag(rX), real(rX));
                A  = 2*(1-rp*rX*wX)/c2;
            else
                A  = (1/X-1.5/X^2)/c2;
            end
        end
    else
        z = ((nu0-nu)*1i+c0)/nuD;
        w = cpf(-imag(z), real(z));
        A = w*rp/nuD;
    end
    I = LM/pi*A/(1-(nuR+NuOptIm*1i)*A);
    
    % Output
    real_part = real(I);
    imag_part = imag(I);
end

