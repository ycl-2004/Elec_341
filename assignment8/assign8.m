%%%%%%%%%%%%%%%%%%%%
% ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

%A = 4 + 10;
%B = 8 + 10;
%C = 8 + 10;
%D = 6 + 10;
%E = 1 + 10;
%F = 7 + 10;
%G = 8 + 10;
%H = 5 + 10;

%SN = 48861785;

A = 60;
B = 10;
C = 40;
D = 40;
E = 50;
F = 50;
G = 70;
H = 40;

SN = 60445574;

s = tf('s');

CF = A*10;
DC = 30/100;
tauf = (B + C)*1e-3;
taui = D*1e-3;
deltaT = 1/CF;
Ntotal = DC+0.5;
% FROM ASSIGN 6

p1 = -1*A;
p2 = -2*B+C*2*i;
p3 = -2*B-C*2*i;
p4 = -3*D;

z1 = -5*E;

p5 = -12*F;

Kp = 3/A;
Ks = 5/B;

Gp = zpk(z1, [p1, p2, p3, p4], 1);
Gp_gain = dcgain(Gp);
Gp_gain_mod = Kp/Gp_gain;

Gp_mod = Gp_gain_mod*Gp;

Hs = zpk([],p5,1);
Hs_gain = dcgain(Hs);
Hs_gain_mod = Ks/Hs_gain;

Hs_mod = Hs_gain_mod*Hs;

Kh = 1/Ks;

GH = Gp_mod*Hs_mod*Kh;

Ntotal = (15*20-42)/65+0.5+DC;
dh = CF/(Ntotal*s+CF);


Q1.G = Gp_mod;
Q1.H = Hs*Kh*Hs_gain_mod*dh;

GH = Q1.G*Q1.H;
ss = stepinfo(GH);

beta = exp(-deltaT/taui);
Q1.Nb =beta/(1-beta)+0.5; 
dh = CF/(Q1.Nb*s+CF);


%----------------- START ASSIGNMENT 8

Q2.D = 1/-z1*(s-z1)/s;

%---------------------------- QUESTION 3 
K_Val =15.3;
Gcl = feedback(K_Val*Q2.D*Q1.G, dh*Q1.H);
%rloc%us(Gcl)
%step(Gcl)
ss = stepinfo(Gcl);
ss
step(Gcl)


Q3.K = K_Val;
Q3.Tr = 0.73; %0.75 2.5

%----------------------------- Question 4
%zeros = zero(Gcl)
%poles = pole(Gcl)

%z = -29.3548;
%p = -34.51;

%Q4.Kp = -1/z/Q3.K;

Q4.Kp  = -1/z1*Q3.K;

%Q4.Ki = (1/p1^2-(z1+p1)/p1*z1)*Q3.K;
Q4.Ki = Q3.K;

%-------------------------

Dp = 1/s;
[Q5.K0,pm,wx,wcp] = margin(Dp*Q1.G*Q1.H);

[K0,pm,Q5.wx,wcp] = margin(Q5.K0*Dp*Q1.G*Q1.H);

%bode(Q5.K0*Dp*Q1.G*Q1.H) % see freqeuncy at -180%

%Q5.K0 = 85;  % 85 IS RIGHT
%Q5.wx = wxo; % 12.5 
%Q5.PM = 29.5;  % Bode plot at || = 1, 

Dnew = (s+Q5.wx)/Q5.wx;
%bode(Q5.K0*Dp*Q1.G*Q1.H*Dnew)  % value -> 25

[K0,Q5.PM,wx,wcp] = margin(Q5.K0*Dp*Q1.G*Q1.H*Dnew);

%-------------------------------


wtest = 1;
pmMax = 1;
while wtest <= 100 % Random number 
    [Gm,Pm,Wcg,Wcp] = margin(Q5.K0 * Dp * Q1.G * Q1.H * (s + wtest) / wtest);
    if pmMax < Pm
        wMax = wtest;
        pmMax = Pm;
    end
    wtest = wtest + 1;
end

Q6.wz = wMax;
Q6.PM = pmMax;

Q6.D = (s-Q6.wz)/(-s*Q6.wz);

%--------------------------------

K_Val = 9.5;
ForwardPath = K_Val*Q1.G*Q6.D;
FeedbackPath = Q1.H;


t = 0:0.001:10;

ClosedLoopSystem = feedback(ForwardPath, FeedbackPath);
%[y, ~] = step(ClosedLoopSystem,t)
%step(ClosedLoopSystem)

%[Q7.K,Pm,Wcg,Wcp] = margin(Q5.K0 * Dp * Q1.G * Q1.H * (s + wMax) / wMax);
Q7.K = K_Val;
Q7.Tr = 1.05;
%-------------------------------

Q8.Kp = Q7.K;
Q8.Ki = Q7.K*dcgain(Q6.D);

%-------------------------------

p = -CF/Q1.Nb;
Dp = -p/(s-p);
[Q9.K0,PM,wx,Wcp] = margin(Dp*Q1.G*Q1.H)
[Gm,PM,Q9.wx,Wcp] = margin(Q9.K0*Dp*Q1.G*Q1.H)

wtest = 1;
pmMax = 1;
while wtest <= 100 % Random number 
    [Gm,Pm,Wcg,Wcp] = margin(Q9.K0 * Q1.G * Q1.H * (s+wtest)/(wtest)*Dp);
    if pmMax < Pm
        wMax = wtest;
        pmMax = Pm;
    end
    wtest = wtest + 1;
end

wMax
Q9.wz = wMax;

%-------------------------------
D10 = Dp*(s+Q9.wz)/Q9.wz;

K_Val = 8.3;

ForwardPath = K_Val*Q1.G*D10;
FeedbackPath = Q1.H;

ClosedLoopSystem = feedback(ForwardPath, FeedbackPath);
%step(ClosedLoopSystem)

Q10.K = K_Val;
Q10.Tr = 0.06;
Q10.Ess = (1-0.642)*100;

%------------------------------------
Q11.Kp = Q10.K;

z = 1/taui;
Q11.Kd = Q10.K*(z-Q9.wz)/(Q9.wz*z)*0.98;

%run a8Submit.p