clear all
global e pi rp sln2 num0 numinf
load('mhtconstants.mat')

nu0     =  0.0;
GamD    = 10.0;
Gam0    =  0.5;
Gam2    =  0.15;
Shift0  =  0.25;
Shift2  =  0.03;
NuOptRe =  0.4;
NuOptIm = -0.05;
x       = -50:0.2:50;
yR      = zeros(1,length(x));
yI      = zeros(1,length(x));

for i=1:length(x)
    [yR(i), yI(i)] = profile(nu0,GamD,Gam0,Gam2,Shift0,Shift2,NuOptRe,NuOptIm,x(i));
end

plot(x,yR)
hold on
plot(x,yI)
