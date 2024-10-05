function [real_part, imag_part] = profilee(nu0, GammaD, Gamma0, Gamma2, Delta0, Delta2, NuOptRe, NuOptIm, nu, varargin)
    % ---------------------------------------- 
    %    PROFILE_mHT: modified Hartman Tran profile
    %    Subroutine to compute the complex normalized spectral shape of an 
    %    isolated line by the mHT model
    %
    %    Standard Input Parameters:
    %    --------------------
    %    nu0       : Unperturbed line position in cm-1.
    %    GammaD    : Doppler HWHM in cm-1.
    %    Gamma0    : Speed-averaged line-width in cm-1.       
    %    Gamma2    : Speed dependence of the line-width in cm-1.
    %    Delta0    : Speed-averaged line-shift in cm-1.
    %    Delta2    : Speed dependence of the line-shift in cm-1.   
    %    NuOptRe   : Real part of the Dicke parameter in cm-1.
    %    NuOptIm   : Imaginary part of the Dicke parameter in cm-1.    
    %    nu        : Current WaveNumber of the Computation in cm-1.
    %
    %    Optional Input Parameters:
    %    --------------------
    %    Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    %    Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    %    alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless (default: 10.0).
    %    dispFlag  : Boolean trigger for including dispersion profile in the output (default: False).
    %
    %    The function has two outputs:
    %    --------------------
    %    (1): Real part of the normalized spectral shape in cm.
    %    (2): Imaginary part of the normalized spectral shape in cm.
    % ---------------------------------------- 
   
    % Definition of standard and optional input parameters
    p   = inputParser;   
    num = @(x) isnumeric(x) && isscalar(x);
    addRequired(p,"nu0",num)
    addRequired(p,"GammaD",num)
    addRequired(p,"Gamma0",num)
    addRequired(p,"Gamma2",num)
    addRequired(p,"Delta0",num)
    addRequired(p,"Delta2",num)
    addRequired(p,"NuOptRe",num)
    addRequired(p,"NuOptIm",num)
    addRequired(p,"nu",num)
    addOptional(p,"Xlm",0.0,num)
    addOptional(p,"Ylm",0.0,num)
    addOptional(p,"alpha",10.0,num)
    parse(p,nu0,GammaD,Gamma0,Gamma2,Delta0,Delta2,NuOptRe,NuOptIm,nu,varargin{:});
    Xlm = p.Results.Xlm;
    Ylm = p.Results.Ylm;
    alpha = p.Results.alpha;

    % Define mathematical, numerical and function handles variables
    pi      = 3.141592653589793;  % Pi number
    rp      = 1.772453850905516;  % Root square of pi
    sln2    = 0.8325546111576977; % Root square of natural logarithm of 2
    num0    = 1.0e-15;            % Numerical zero
    numinf  = 4000;               % Numerical infinity
    cpf     = @cpf_accurate;   % cpf choice
    ternary = @(varargin) varargin{end - varargin{1}}; % compressed if-else statement

    % The function itself
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

    real_part = real(I);
    imag_part = imag(I);
end
