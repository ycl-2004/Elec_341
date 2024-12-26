%%%%%%%%%%%%%%%%%%%%
% Assignment 4 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;


%--------------------------------------
% Question 1 

syms U1 U2 Y1 Y2 s;

G1 = 1/(s+ (4+10));
G2 = 10/(s+ (8+10));
G3 = 10/(s+ (8+10));
G4 = 10/(s+ (6+10));
G5 = 1/(s+ (1+10));

H1 = 100/(s+ (7+10));
H2 = 500/(s+ (8+10));

U2 = 0;

eq1 = Y1/G2 == (U1-Y2*H2)*G1;
eq3 = Y2 == G5*((Y1*H1-U2)*G4+(U1-Y2*H2)*G1*G3);
%eq2 = A*G2 == Y1;

[y1,y2] = solve([eq1,eq3],[Y1,Y2]);

h11 = simplify(y1/U1);
[num, den] = numden(h11);
h11_final = tf(sym2poly(num), sym2poly(den));

h21 = simplify(y2/U1);
[num, den] = numden(h21);
h21_final = tf(sym2poly(num), sym2poly(den));

clear U2 Y1 Y2 y1 y2; 
syms U2 Y1 Y2; 

U1 = 0;

eq1 = Y1/G2 == (U1-Y2*H2)*G1;
eq3 = Y2 == G5*((Y1*H1-U2)*G4+(U1-Y2*H2)*G1*G3);

[y1,y2] = solve([eq1,eq3],[Y1,Y2]);

h12 = simplify(y1/U2);
[num, den] = numden(h12);
h12_final = tf(sym2poly(num), sym2poly(den));

h22 = simplify(y2/U2);
[num, den] = numden(h22);
h22_final = tf(sym2poly(num), sym2poly(den));


Q1.G11 = h11_final;
Q1.G12 = h12_final;
Q1.G21 = h21_final;
Q1.G22 = h22_final;

% Question 1 End
%--------------------------------------
% Question 2 


clear U2 Y1 Y2 y1 y2; 
syms U Y1 Y2; 

eq1 = Y1/G2 == (U-Y2*H2)*G1;
eq2 = Y2 == G5*((Y1*H1-U)*G4+(U-Y2*H2)*G1*G3);

[y1,y2] = solve([eq1,eq2],[Y1,Y2]);

h1 = simplify(y1/U);
[num, den] = numden(h1);
h1_final = tf(sym2poly(num), sym2poly(den));

h2 = simplify(y2/U);
[num, den] = numden(h2);
h2_final = tf(sym2poly(num), sym2poly(den));

Q2.G1 = h1_final;
Q2.G2 = h2_final;


% Question 2 End
%--------------------------------------
% Question 3

%Assignment 3 Parameter 
clear s;
s = tf("s");

Rw = (4+10)/3;
Lw = (8+10)*1e-3;
Jr = (8+10)/10*1e-3;
Br = ((6+10)+(1+10))*1e-3;
Ke = (8+10)*50*1e-3;
Km = (8+10)*50*1e-3;
%------------------------

Jf = (4+10)/25*1e-3;
Bf = (8+10)/20*1e-3;
Jg = (8+10)*50*1e-6;
Bg = (6+10)/30*1e-3;
n = (1+10);

Jtotal = Jr + Jg*n^2 + Jf*n^2;
Btotal = Br + Bg*n^2 + Bf*n^2;

Q3.Ye = minreal(1/(s*Lw + Rw));
Q3.Ym = minreal(1/(s*Jtotal + Btotal));

Q4.Gi = feedback(Q3.Ye,Km*Q3.Ym*Km); %Gi = iw/Vin
Q4.Gw = feedback(Q3.Ye*Km*Q3.Ym,Km)*n; % Gw = Wr/Vin
run a4Submit.p
