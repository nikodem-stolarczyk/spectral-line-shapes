function res = cpf_fast(x, y)
    % CPF_FAST: Computes the complex probability function using Humlicek's 
    % algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) 
    % and using a rational series with 24 terms in other subregions. (Source: jstor.org/stable/2158232)
    %
    % Input/Output Parameters of Routine
    % x : Real part of input complex parameter
    % y : Imaginary part of input complex parameter
    %
    % The function has one output:
    % (1): Complex probability function
    global rp
    hum1_threshold = 15.0; % the border of the first Humlicek's region

    if abs(x) + y > hum1_threshold
        t   = y - x * 1i;
        res =  t / (0.5 + t.^2) / rp;
    else
        z   = -y + x*1i;
        L   = 4.119534287814236;
        Z   = (L + z) / (L - z);
        a   = [ -1.513747622620502E-10,  4.904820407381768E-09,  1.331045329581992E-09, -3.008282344381996E-08, ...
                -1.912225887484805E-08,  1.873834346505099E-07,  2.568264135399530E-07, -1.085647579417637E-06, ...
                -3.038893184366094E-06,  4.139461724429617E-06,  3.047106608295325E-05,  2.433141546207148E-05, ...
                -2.074843151143828E-04, -7.816642995626165E-04, -4.936426901286291E-04,  6.215006362949147E-03, ...
                 3.372336685531603E-02,  1.083872348456673E-01,  2.654963959880772E-01,  5.361139535729116E-01, ...
                 9.257087138588670E-01,  1.394819673379119E+00,  1.856286499205540E+00,  2.197858936531542E+00];
        res = 2*sum(a.*Z.^(23:-1:0))/(L-z)^2 + 1/(L-z)/rp;
    end

end
