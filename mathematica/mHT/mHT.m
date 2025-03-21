(* ::Package:: *)

BeginPackage["mHT`"];


CPFAccurate::usage      = "Computes the complex probability function using accurate method.";
CPFFast::usage          = "Computes the complex probability function using fast method.";
mHTProfile::usage       = "Computes modified Hartmann Tran profile function.";
mHTProfileVector::usage = "Computer modified Hartmann Tran profile function with frequency as table input.";


Begin["`Private`"];


(* ----------------------------------------
    Accurate CPF algorithm
    Computes the complex probability function using a rational series
    with 42 terms. It is assumed that Im(z) > 0 or Im(z) = 0. (reference: jstor.org/stable/2158232)
    A series was simplified to 37 terms introducing less than 10^(-17)
    deviations on mHT profile.
    
    To decrease code execution time, following 
    numeric values were introduced explicitly:
     * Sqrt[1/Pi] = 0.5641895835477563`
    
    Standard Input Parameters:
    -------------------- 
    x  : Real part of input complex parameter
    y  : Imaginary part of input complex parameter
    
    The function has one output:
    --------------------
    (1): Complex probability function
   ---------------------------------------- *) 
CPFAccurate[{x_, y_}] := 
Block[{
  z   = -y + I x,
  L   = 5.449631621480024` (* The pre-calculated Weideman constant for N = 42 *),
  a   = { -3.129493160727961`*^-14, -1.188364999909099`*^-14, +1.951777029849348`*^-13, +1.790586243645278`*^-13,
          -1.184560208678836`*^-12, -2.069163661083667`*^-12, +6.430136110306704`*^-12, +2.063579921011804`*^-11,
          -2.392389527320517`*^-11, -1.799169607159564`*^-10, -6.353807951660892`*^-11, +1.282896083944607`*^-09,
          +2.636162411919059`*^-09, -5.468780625369738`*^-09, -3.294773119114329`*^-08, -2.752070035718561`*^-08,
          +2.206733163926054`*^-07, +8.511689670641750`*^-07, +4.936972061734341`*^-07, -6.617492208403963`*^-06,
          -2.914574364851397`*^-05, -4.816473680511106`*^-05, +1.044072210002090`*^-04, +1.070131083157417`*^-03,
          +4.631075611097791`*^-03, +1.480296368764821`*^-02, +3.922970169744468`*^-02, +9.038744880336540`*^-02,
          +1.857036333535562`*^-01, +3.455278077566057`*^-01, +5.882708203344523`*^-01, +9.230959991941070`*^-01,
          +1.342044484596932`*^+00, +1.814714451499866`*^+00, +2.288734169675538`*^+00, +2.697763665856064`*^+00,
          +2.975931371735470`*^+00} (* The pre-calculated table of FFT constant terms (truncated from the begining) *)},
  (2 a . ((L+z)/(L-z))^Range[36,0,-1] / (L-z) + 0.5641895835477563`) / (L-z)] /; 
NumericQ[x] && NumericQ[y]


(* ----------------------------------------
    Fast CPF algorithm
    Computes the complex probability function using Humlicek's 
    algorithm in its first subregion (Source: 10.1016/0022-4073(82)90078-4) 
    and using a rational series with 24 terms in other subregions. (reference: jstor.org/stable/2158232)
    
    To decrease code execution time, following 
    numeric values were introduced explicitly:
     * Sqrt[1/Pi] = 0.5641895835477563`
     
    Standard Input Parameters:
    -------------------- 
    x  : Real part of input complex parameter
    y  : Imaginary part of input complex parameter
    
    The function has one output:
    --------------------
    (1): Complex probability function
   ---------------------------------------- *) 
CPFFast[{x_, y_}]:=
Block[{
  z   = -y + I x, 
  L   = 4.119534287814236` (* The pre-calculated Weideman constant for N = 24 *),
  a   = { -1.513747622620502`*^-10,  4.904820407381768`*^-09,  1.331045329581992`*^-09, -3.008282344381996`*^-08,
          -1.912225887484805`*^-08,  1.873834346505099`*^-07,  2.568264135399530`*^-07, -1.085647579417637`*^-06,
          -3.038893184366094`*^-06,  4.139461724429617`*^-06,  3.047106608295325`*^-05,  2.433141546207148`*^-05,
          -2.074843151143828`*^-04, -7.816642995626165`*^-04, -4.936426901286291`*^-04,  6.215006362949147`*^-03,
           3.372336685531603`*^-02,  1.083872348456673`*^-01,  2.654963959880772`*^-01,  5.361139535729116`*^-01,
           9.257087138588670`*^-01,  1.394819673379119`*^+00,  1.856286499205540`*^+00,  
           2.197858936531542`*^+00} (* The pre-calculated table of FFT constant terms *)},
  If[Abs[x] + y > 15.0 (* The border of the first Humlicek's region *), 
    -0.5641895835477563` z / (z^2+0.5),
    (2 a . ((L+z)/(L-z))^Range[23,0,-1] / (L-z) + 0.5641895835477563`) / (L-z)]]/;
NumericQ[x] && NumericQ[y] 


(* ----------------------------------------
    Beta-Correction function
    Subroutine to compute beta-correction used for hard-collision based line-shape profiles
    To correct NuOptRe value in the profile . Applicable up to alpha = 5.0, for higher 
    alpha values - correction neglected. 
    Source: 10.1016/j.jqsrt.2019.106784
    
    Standard Input Parameters:
    --------------------
    GammaD  : Doppler broadening in cm-1. 
    NuOptRe : Real part of the Dicke parameter in cm-1.
    alpha   : Mass ratio in the molecule, applicable up to alpha = 5.0, dimensionless.
   
    The function has one output:
    --------------------
    (1)     : Value of the beta correction, dimensionless. 
   ---------------------------------------- *)
\[Beta][{GammaD_, NuOptRe_, alpha_}] :=
If[alpha < 5.0 (* Maximum applicable alpha value *),
  (* a Tanh[b Log10[NuOptRe/GammaD] + c] + d *)
  (* a *)(0.0534` + 0.1585` Exp[-0.4510` alpha])\
  *Tanh[\
  (* b *)((((0.0050` alpha) + 0.0056`) alpha - 0.1258`) alpha + 1.9595`)\
  Log10[NuOptRe/GammaD] +\
  (* c *)((((0.0003` alpha) - 0.0125`) alpha + 0.0672`) alpha - 0.0546`)\
  ]+\
  (* d *)(0.9466` - 0.1585` Exp[-0.4510` alpha]),
  1.0] /;
NumericQ[GammaD] && NumericQ[NuOptRe] && NumericQ[alpha] 


(* Chosing the CPF used for the profile calculation *)
CPFChoice = 1;
CPF      := ({CPFAccurate,CPFFast})[[CPFChoice]];


(* ----------------------------------------
    Modified Hartmann Tran profile function
	Subroutine to compute the complex normalized spectral-line shape using mHT model.
    
    To decrease code execution time, following 
    numeric values were introduced explicitly:
	 * 1/Pi           = 0.3183098861837907`
	 * Sqrt[Pi]       = 1.772453850905516`
	 * 2 Sqrt[Pi]     = 3.5449077018110318`
	 * 1/Sqrt[Log[2]] = 1.2011224087864498`
    
    Standard Input Parameters:
    --------------------
    nu0       : Unperturbed line position in cm-1.
    GammaD    : Doppler broadening in cm-1.
    Gamma0    : Speed-averaged line-width in cm-1.       
    Gamma2    : Quadratic speed dependence parameter of the line-width in cm-1.
    Delta0    : Speed-averaged line-shift in cm-1.
    Delta2    : Quadratic speed dependence parameter of the line-shift in cm-1.   
    NuOptRe   : Real part of the complex Dicke parameter in cm-1.
    NuOptIm   : Imaginary part of the complex Dicke parameter in cm-1.    
    nu        : Current wavenumber of the computation in cm-1.
    
    Optional Input Parameters:
    --------------------
    Ylm       : Imaginary part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    Xlm       : Real part of the 1st order (Rosenkranz) line mixing coefficients, dimensionless (default: 0.0).
    alpha     : Mass ratio in the molecule for calculating beta-correction, applicable up to alpha=5, dimensionless (default: 10.0).
    disp      : Boolean trigger for including dispersion profile in the output (default: False).
    
    The function has one outputs:
    --------------------
    (1)       : Real or imaginary (depending on disp value) part of the normalized spectral shape in cm.
   ---------------------------------------- *)
mHTProfile[{nu0_,GammaD_,Gamma0_,Gamma2_,Delta0_,Delta2_,NuOptRe_,NuOptIm_,nu_,Ylm_:0,Xlm_:0,alpha_:10,disp_:False}] :=
Block[{nuD = 1.2011224087864498` GammaD,
    nuR,c0,c2,LM,X,Y,csqY,z1,z2,w1,w2,rX,wX,z,w,A,ImHT},
    nuR = NuOptRe \[Beta][{GammaD,NuOptRe,alpha}];
    c2  = Gamma2 + I Delta2;
    c0  = Gamma0 + I Delta0 - 1.5 c2 + nuR + I NuOptIm;
    LM  = 1 + Xlm - I Ylm;
    If[Abs[c2] > (1.0`*^-9) (* Limit where speed dependence impact is lower than numerical noise level *), 
        X    = (I (nu0 - nu) + c0) / c2;
        Y    = 0.25 (nuD / c2)^2;
        csqY = 0.50 nuD (Gamma2 - I Delta2) / (Gamma2^2+ Delta2^2);
        If[Abs[Y]>Abs[X * (1.0`*^-15) (* Numerical zero *)] ,
            z2 = Sqrt[X + Y] + csqY;
            z1 = If[Abs[X] > Abs[Y * (3.`*^-8)], z2 - 2 csqY, (I (nu0 - nu) + c0) / nuD ];
            w1 = CPF[ReIm[I z1]];
            w2 = CPF[ReIm[I z2]];
            A  = 1.772453850905516` / nuD (w1 - w2);,
            If[Abs[Sqrt[X]] < (4.`*^3) (* Numerical infinity *),
                rX = Sqrt[X];
                A  = (2 - 3.5449077018110318` rX CPF[ReIm[I rX]]) / c2;,
                A  = (1 - 1.5 / X) / X / c2;
            ];
        ],
        z = (I (nu0 - nu) + c0) / nuD;
        w = CPF[ReIm[I z]];
        A = w 1.772453850905516` / nuD;
    ];
    ImHT = 0.3183098861837907` LM A / (1 - A (nuR + I NuOptIm));
    If[Not[disp],Re[ImHT],-Im[ImHT]]
]/;
NumericQ[nu0] && NumericQ[GammaD] && NumericQ[Gamma0] && NumericQ[Gamma2] && NumericQ[Delta0] && NumericQ[Delta2] &&\
NumericQ[NuOptRe] && NumericQ[NuOptIm] && NumericQ[Ylm] && NumericQ[Xlm] && NumericQ[alpha]


(* Accurate CPF algorithm rewritted to allow vector input *)
CPFAccurateVector[xytab_] :=
Block[{
  z   = -xytab[[All,2]]+ I xytab[[All,1]],
  L   = 5.449631621480024` (* The pre-calculated Weideman constant for N = 42 *),
  a   = { -3.129493160727961`*^-14, -1.188364999909099`*^-14, +1.951777029849348`*^-13, +1.790586243645278`*^-13,
          -1.184560208678836`*^-12, -2.069163661083667`*^-12, +6.430136110306704`*^-12, +2.063579921011804`*^-11,
          -2.392389527320517`*^-11, -1.799169607159564`*^-10, -6.353807951660892`*^-11, +1.282896083944607`*^-09,
          +2.636162411919059`*^-09, -5.468780625369738`*^-09, -3.294773119114329`*^-08, -2.752070035718561`*^-08,
          +2.206733163926054`*^-07, +8.511689670641750`*^-07, +4.936972061734341`*^-07, -6.617492208403963`*^-06,
          -2.914574364851397`*^-05, -4.816473680511106`*^-05, +1.044072210002090`*^-04, +1.070131083157417`*^-03,
          +4.631075611097791`*^-03, +1.480296368764821`*^-02, +3.922970169744468`*^-02, +9.038744880336540`*^-02,
          +1.857036333535562`*^-01, +3.455278077566057`*^-01, +5.882708203344523`*^-01, +9.230959991941070`*^-01,
          +1.342044484596932`*^+00, +1.814714451499866`*^+00, +2.288734169675538`*^+00, +2.697763665856064`*^+00,
          +2.975931371735470`*^+00} (* The pre-calculated table of FFT constant terms (truncated from the begining) *)},
  Table[(2 a . ((L+z[[j]])/(L-z[[j]]))^Range[36,0,-1] / (L-z[[j]]) + 0.5641895835477563`) / (L-z[[j]]),{j,1,Length[z]}]]


(* Fast CPF algorithm rewritted to allow vector input *)
CPFFastVector[xytab_]:=
Block[{
  z   = -y + I x, 
  L   = 4.119534287814236` (* The pre-calculated Weideman constant for N = 24 *),
  a   = { -1.513747622620502`*^-10,  4.904820407381768`*^-09,  1.331045329581992`*^-09, -3.008282344381996`*^-08,
          -1.912225887484805`*^-08,  1.873834346505099`*^-07,  2.568264135399530`*^-07, -1.085647579417637`*^-06,
          -3.038893184366094`*^-06,  4.139461724429617`*^-06,  3.047106608295325`*^-05,  2.433141546207148`*^-05,
          -2.074843151143828`*^-04, -7.816642995626165`*^-04, -4.936426901286291`*^-04,  6.215006362949147`*^-03,
           3.372336685531603`*^-02,  1.083872348456673`*^-01,  2.654963959880772`*^-01,  5.361139535729116`*^-01,
           9.257087138588670`*^-01,  1.394819673379119`*^+00,  1.856286499205540`*^+00,  
           2.197858936531542`*^+00} (* The pre-calculated table of FFT constant terms *)},
  If[Abs[x] + y > 15.0 (* The border of the first Humlicek's region *), 
    -0.5641895835477563` z / (z^2+0.5),
    Table[(2 a . ((L+z[[j]])/(L-z[[j]]))^Range[23,0,-1] / (L-z[[j]]) + 0.5641895835477563`) / (L-z[[j]]),{j,1,Length[z]}]]]
  
  
(* Chosing the CPF used for the profile calculation - consistent with CPF choice *)
CPFVector := ({CPFAccurateVector,CPFFastVector})[[CPFChoice]];


mHTProfileVector[{nu0_,GammaD_,Gamma0_,Gamma2_,Delta0_,Delta2_,NuOptRe_,NuOptIm_,nu_,Ylm_:0,Xlm_:0,alpha_:10,disp_:False}] :=
Block[{nuD = 1.2011224087864498` GammaD,
    nuR,c0,c2,LM,X,Y,csqY,z1,z2,w1,w2,rX,wX,z,w,A,ImHT},
    nuR = NuOptRe \[Beta][{GammaD,NuOptRe,alpha}];
    c2  = Gamma2 + I Delta2;
    c0  = Gamma0 + I Delta0 - 1.5 c2 + nuR + I NuOptIm;
    LM  = 1 + Xlm - I Ylm;
    If[Abs[c2] > (1.0`*^-9) (* Limit where speed dependence impact is lower than numerical noise level *), 
        X    = (I (nu0 - nu) + c0) / c2;
        Y    = 0.25 (nuD / c2)^2;
        csqY = 0.50 nuD (Gamma2 - I Delta2) / (Gamma2^2+ Delta2^2);
        A    = Table[If[Abs[Y]>Abs[X[[j]] * (1.0`*^-15) (* Numerical zero *)] ,
            z2 = Sqrt[X[[j]] + Y] + csqY;
            z1 = If[Abs[X[[j]]] > Abs[Y * (3.`*^-8)], z2 - 2 csqY, (I (nu0 - nu[[j]]) + c0) / nuD ];
            w1 = CPF[ReIm[I z1]];
            w2 = CPF[ReIm[I z2]];
            1.772453850905516` / nuD (w1 - w2),
            If[Abs[Sqrt[X[[j]]]] < (4.`*^3) (* Numerical infinity *),
                rX = Sqrt[X[[j]]];
                (2 - 3.5449077018110318` rX CPFv[ReIm[I rX]]) / c2,
                (1 - 1.5 / X[[j]]) / X[[j]] / c2
            ];
        ],{j,1,Length[X]}];,
        z = (I (nu0 - nu) + c0) / nuD;
        w = CPFVector[ReIm[I z]];
        A = w 1.772453850905516` / nuD;
    ];
    ImHT = 0.3183098861837907` LM A / (1 - A (nuR + I NuOptIm));
    If[Not[disp],Re[ImHT],-Im[ImHT]]
]/;
NumericQ[nu0] && NumericQ[GammaD] && NumericQ[Gamma0] && NumericQ[Gamma2] && NumericQ[Delta0] && NumericQ[Delta2] &&\
NumericQ[NuOptRe] && NumericQ[NuOptIm] && NumericQ[Ylm] && NumericQ[Xlm] && NumericQ[alpha]


End[ ];


EndPackage[ ];
