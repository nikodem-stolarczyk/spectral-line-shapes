function res = cpf_fast_vector(x, y)
    % ---------------------------------------- 
    %    Computes the complex probability function using Humlicek's 
    %    algorithm in its first subregion (reference: 10.1016/0022-4073(82)90078-4) 
    %    and using a rational series with 24 terms in other subregions. (reference: jstor.org/stable/2158232)
    %
    %    To decrease code execution time, following 
    %    numeric values were introduced explicitly:
    %    * Sqrt[1/Pi] = 0.5641895835477563
    %
    %    Standard Input Parameters:
    %    --------------------
    %    x   : Real part of input complex parameter.
    %    y   : Imaginary part of input complex parameter.
    %
    %    The function has one output:
    %    --------------------
    %    (1) : Complex probability function.
    % ---------------------------------------- 

    if abs(x) + y > 15.0 % The border of the first Humlicek's region
        t   = y - x * 1i;
        res =  0.5641895835477563 * t ./ (0.5 + t.^2) ;
    else
        z   = -y + x*1i;
        L   = 4.119534287814236; % The pre-calculated Weideman constant for N = 24
        Z   = (L + z) ./ (L - z);
        a   = [ -1.513747622620502E-10,  4.904820407381768E-09,  1.331045329581992E-09, -3.008282344381996E-08, ...
                -1.912225887484805E-08,  1.873834346505099E-07,  2.568264135399530E-07, -1.085647579417637E-06, ...
                -3.038893184366094E-06,  4.139461724429617E-06,  3.047106608295325E-05,  2.433141546207148E-05, ...
                -2.074843151143828E-04, -7.816642995626165E-04, -4.936426901286291E-04,  6.215006362949147E-03, ...
                 3.372336685531603E-02,  1.083872348456673E-01,  2.654963959880772E-01,  5.361139535729116E-01, ...
                 9.257087138588670E-01,  1.394819673379119E+00,  1.856286499205540E+00,  2.197858936531542E+00]; % The pre-calculated table of FFT constant terms
        res = (2*sum(a.*(Z.'.^(23:-1:0)),2)./(L-z).'+ 0.5641895835477563)./(L-z).';
    end
end
