function w42 = cpf_accurate(x, y)
    % ---------------------------------------- 
    %    CPF_ACCURATE: Computes the complex probability function using a rational series 
    %    with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0. (Source: jstor.org/stable/2158232)
    %    A series was simplified to 37 terms introducing less than 10^(-17)
    %    deviations on mHT profile.
    %
    %    Standard Input Parameters:
    %    --------------------
    %    x : Real part of input complex parameter
    %    y : Imaginary part of input complex parameter
    %
    %    The function has one output:
    %    --------------------
    %    (1): Complex probability function
    % ---------------------------------------- 

    rp  = 1.772453850905516; % root square of pi 
    L   = 5.449631621480024; % the pre-calculated Weideman constant for N = 42
    z   = -y + x*1i;
    Z   = (L + z) / (L - z);
    a   = [ -3.129493160727961E-14, -1.188364999909099E-14,  1.951777029849348E-13,  1.790586243645278E-13, ...
            -1.184560208678836E-12, -2.069163661083667E-12,  6.430136110306704E-12,  2.063579921011804E-11, ...
            -2.392389527320517E-11, -1.799169607159564E-10, -6.353807951660892E-11,  1.282896083944607E-09, ...
             2.636162411919059E-09, -5.468780625369738E-09, -3.294773119114329E-08, -2.752070035718561E-08, ...
             2.206733163926054E-07,  8.511689670641750E-07,  4.936972061734341E-07, -6.617492208403963E-06, ...
            -2.914574364851397E-05, -4.816473680511106E-05,  1.044072210002090E-04,  1.070131083157417E-03, ...
             4.631075611097791E-03,  1.480296368764821E-02,  3.922970169744468E-02,  9.038744880336540E-02, ...
             1.857036333535562E-01,  3.455278077566057E-01,  5.882708203344523E-01,  9.230959991941070E-01, ...
             1.342044484596932E+00,  1.814714451499866E+00,  2.288734169675538E+00,  2.697763665856064E+00, ...
             2.975931371735470E+00]; % the pre-calculated table of FFT constant terms (truncated from the beginning)

    w42 = 2*sum(a.*Z.^(36:-1:0))./(L-z).^2 + 1./(L-z)./rp;
end
