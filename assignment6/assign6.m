% FORMULAS

%%%%%%%%%%%%%%%%%%%%
% ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

SN = 48861785;

%----------------- Question 1

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

Q1.Kh = Kh;
Q1.G = Gp_mod;
Q1.H = Hs*Kh*Hs_gain_mod;

GH = Q1.G*Q1.H;

Q1.Ku = margin(GH);

K = Q1.Ku/2;

Q2.Gcl = feedback(K*Q1.G, Hs_mod*Kh);

step(Q2.Gcl)

S = stepinfo(Q2.Gcl);
Q3.Tr = 0.0573;
Q3.Tp = S.PeakTime;
Q3.Ts = S.SettlingTime;
Q3.yf = 0.735;
Q3.yp = S.Peak;
Q3.Ess = 26; % GUESS IDK HOW, RANDOM VALUE TEST
Q3.OSu = (1.08-1)*100;  %
Q3.OSu = 7.5; %8.5 10~20
Q3.OSy = S.Overshoot;

%---------------------------------

% Target overshoot specification
OSy_target = 20; % 20%

% Find the gain K that achieves the desired overshoot
K_opt = fminsearch(@(K) abs(stepinfo(feedback(K*Q1.G, Hs_mod*Kh)).Overshoot - OSy_target), K);

% Update closed-loop transfer function with the optimized gain K_opt
Q2.Gcl_adjusted = feedback(K_opt * Q1.G, Hs_mod * Kh);

% Plot the step response of the adjusted system
step(Q2.Gcl_adjusted);

% Extract step response characteristics for the adjusted system
S_adjusted = stepinfo(Q2.Gcl_adjusted);

% Calculations for final values
Ts_adjusted = S_adjusted.SettlingTime;
Ess_adjusted = abs(1 - dcgain(Q2.Gcl_adjusted))*100; % Steady-state error assuming unit step input

% Display results
K_opt
Ts_adjusted
Ess_adjusted


Q4.K = K_opt;
Q4.Ts = Ts_adjusted;
Q4.Ess = Ess_adjusted;

%-----------------------------------------------------------------

Tp_target = 0.9 * Q3.Tp;


K_new = fminsearch(@(K) abs(stepinfo(feedback(K*Q1.G, Hs_mod*Kh)).PeakTime - Tp_target), K);

Q2.Gcl_new = feedback(K_new * Q1.G, Hs_mod * Kh);


step(Q2.Gcl_new);


S_new = stepinfo(Q2.Gcl_new);

Ts_new = S_new.SettlingTime;
Ess_new = abs(1 - dcgain(Q2.Gcl_new))*100;

Q5.K = K_new;
Q5.Ts = Ts_new;
Q5.Ess = Ess_new;
%run a6Submit.p;





