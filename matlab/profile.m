function [real_part, imag_part] = profile(nu0, GamD, Gam0, Gam2, Shift0, Shift2, NuOptRe, NuOptIm, nu, Sw, Ylm, Xlm, alpha)
    % PROFILE_mHT: modified Hartman Tran profile
    % Subroutine to compute the complex normalized spectral shape of an 
    % isolated line by the mHT model
    %
    % Input/Output Parameters of Routine (Arguments or Common)
    % nu0       : Unperturbed line position in cm-1 (Input).
    % GamD      : Doppler HWHM in cm-1 (Input)
    % Gam0      : Speed-averaged line-width in cm-1 (Input).       
    % Gam2      : Speed dependence of the line-width in cm-1 (Input).
    % Shift0    : Speed-averaged line-shift in cm-1 (Input).
    % Shift2    : Speed dependence of the line-shift in cm-1 (Input)   
    % NuOptRe   : Real part of the Dicke parameter in cm-1 (Input).
    % NuOptIm   : Imaginary part of the Dicke parameter in cm-1 (Input).    
    % nu        : Current WaveNumber of the Computation in cm-1 (Input).
    % Sw        : Statistical weight 
    % Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients in cm-1 (Input)
    % Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients in cm-1 (Input)
    % alpha     : Mass ratio in the molecule for calculating beta-correction. Applicable up to alpha=5.
    %
    % The function has two outputs:
    % (1): Real part of the normalized spectral shape (cm)
    % (2): Imaginary part of the normalized spectral shape (cm)
    global e pi rp sln2 num0 numinf
    if nargin < 14
        alpha = 10.0;
        if nargin < 13
            Xlm = 0.0;
            if nargin < 12
                Ylm = 0.0;
                if nargin < 11
                    Sw = 1.0;
                end
            end
        end
    end
    ternary = @(varargin) varargin{end - varargin{1}};
    cpf     = @cpf_accurate;
    

    nuD = GamD / sln2 ;
    nuR = NuOptRe*beta(GamD, NuOptRe, alpha);
    c2  = Gam2+Shift2*1i;
    c0  = Gam0+Shift0*1i-1.5*c2+nuR+NuOptIm*1i;
    LM  = 1+Xlm+Ylm*1i;

    if abs(c2) ~= num0
        X    = ((nu0-nu)*1i+c0)/c2;
        Y    = 0.25*(nuD/c2)^2;
        csqY = 0.50*nuD*(Gam2-Shift2*1i)/(Gam2^2+Shift2^2);
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
    I = Sw*LM/pi*A/(1-(nuR+NuOptIm*1i)*A);

    real_part = real(I);
    imag_part = imag(I);
end


