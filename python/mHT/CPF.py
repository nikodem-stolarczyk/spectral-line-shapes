inverse_sqrt_pi = 0.5641895835477563

def cpf_accurate(x,y):
    """    
    # ----------------------------------------
    #      "CPF_ACCURATE": Accurate CPF algorithm
    #      Computes the complex probability function using a rational series 
    #      with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0. (reference: jstor.org/stable/2158232)
    #      A series was simplified to 37 terms introducing less than 10^(-17)
    #      deviations on mHT profile.
    #
    #      Standard Input Parameters:
    #      -------------------- 
    #      x  : Real part of input complex parameter
    #      y  : Imaginary part of input complex parameter
    #
    #      The function has one output:
    #      --------------------
    #      (1): Complex probability function
    # ----------------------------------------
    """
    z = -y + x*1j
    L = 5.449631621480024 # The pre-calculated Weideman constant for N = 42
    Z = (L+z)/(L-z)
    a = [ -3.129493160727961E-14, -1.188364999909099E-14,  1.951777029849348E-13,  1.790586243645278E-13, 
          -1.184560208678836E-12, -2.069163661083667E-12,  6.430136110306704E-12,  2.063579921011804E-11, 
          -2.392389527320517E-11, -1.799169607159564E-10, -6.353807951660892E-11,  1.282896083944607E-09,  
           2.636162411919059E-09, -5.468780625369738E-09, -3.294773119114329E-08, -2.752070035718561E-08,  
           2.206733163926054E-07,  8.511689670641750E-07,  4.936972061734341E-07, -6.617492208403963E-06, 
          -2.914574364851397E-05, -4.816473680511106E-05,  1.044072210002090E-04,  1.070131083157417E-03, 
           4.631075611097791E-03,  1.480296368764821E-02,  3.922970169744468E-02,  9.038744880336540E-02, 
           1.857036333535562E-01,  3.455278077566057E-01,  5.882708203344523E-01,  9.230959991941070E-01, 
           1.342044484596932E+00,  1.814714451499866E+00,  2.288734169675538E+00,  2.697763665856064E+00, 
           2.975931371735470E+00] # The pre-calculated table of FFT constant terms (truncated from the begining)
    p= a[36] + (a[35] + (a[34] + (a[33] + (a[32] + (a[31] + (a[30] + (a[29] + (a[28] + (a[27] + (a[26] + (a[25] + (a[24] + (a[23] + (a[22] + (a[21] + (a[20] + (a[19] + (a[18] + (a[17] + (a[16] + (a[15] + (a[14] + (a[13] +(a[12] +  (a[11] + (a[10] + (a[9] + (a[8] + (a[7] + (a[6] + (a[5] + (a[4] + (a[3] + (a[2] + (a[1] + a[0] * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z)* Z
    return (2*p/(L-z) + inverse_sqrt_pi)/(L-z)

def cpf_fast(x,y):
    """    
    # ----------------------------------------
    #      "CPF_FAST": Fast CPF algorithm
    #      Computes the complex probability function using Humlicek's 
    #      algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) 
    #      and using a rational series with 24 terms in other subregions. (reference: jstor.org/stable/2158232)
    #
    #      Standard Input Parameters:
    #      -------------------- 
    #      x  : Real part of input complex parameter
    #      y  : Imaginary part of input complex parameter
    #
    #      The function has one output:
    #      --------------------
    #      (1): Complex probability function
    # ----------------------------------------
    """
    hum1_threshold = 15.0 # the border of the first Humlicek's region
    if abs(x)+y>hum1_threshold:
        t = y-x*1j
        return inverse_sqrt_pi*t/(0.5+pow(t,2))
    else: 
        z = -y + x*1j
        L = 4.119534287814236 # the pre-calculated Weideman constant for N = 24
        Z = (L+z)/(L-z)
        a = [ -1.513747622620502E-10,  4.904820407381768E-09,  1.331045329581992E-09, -3.008282344381996E-08,
              -1.912225887484805E-08,  1.873834346505099E-07,  2.568264135399530E-07, -1.085647579417637E-06,
              -3.038893184366094E-06,  4.139461724429617E-06,  3.047106608295325E-05,  2.433141546207148E-05,
              -2.074843151143828E-04, -7.816642995626165E-04, -4.936426901286291E-04,  6.215006362949147E-03,
               3.372336685531603E-02,  1.083872348456673E-01,  2.654963959880772E-01,  5.361139535729116E-01,
               9.257087138588670E-01,  1.394819673379119E+00,  1.856286499205540E+00,  2.197858936531542E+00] 
            # the pre-calculated table of FFT constant terms
        p = a[23] + (a[22] + (a[21] + (a[20] + (a[19] + (a[18] + (a[17] + (a[16] + (a[15] + (a[14] + (a[13] +(a[12] +  (a[11] + (a[10] + (a[9] + (a[8] + (a[7] + (a[6] + (a[5] + (a[4] + (a[3] + (a[2] + (a[1] + a[0] * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z) * Z
        return (2*p/(L-z) + inverse_sqrt_pi)/(L-z)
