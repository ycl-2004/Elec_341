%%%%%%%%%%%%%%%%%%%%
% Assignment 4 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;

n = 2;

M1 = (4+10)*0.001;
M2 = (8+10)*0.001/n^2;
M3 = (8+10)*3*0.001/n^2;

B1 = (6+10)/100;
B23 = (1+10)/500/n^2;

K2 = (7+10)/n^2;
K23 = (8+10)/n^2;



C1 = M1;
C2 = M2;
C3 = M3;

R1 = 1/(B1);
R23 = 1/B23;

L2 = 1/K2;
L23 = 1/K23;

F0 = 1;

s = tf('s');


%CORRECT PART HERE ------------------------------------- 

Y = [(C1+C2)*s+1/R1+1/(L2*s)+1/R23+1/(L23*s) -1/R23-1/(L23*s)
    -1/R23-1/(L23*s) C3*s+1/R23+1/(L23*s)];
I = [F0;0];
V = Y\I;
V = V/s;
Q1.G = V(1);
Q1.G

figure;
step(Q1.G);
title('Step Response of Q1.G');
grid on;

% If you want to get the step response data
[y, t] = step(Q1.G);
S = stepinfo(Q1.G);

S

Q2.zeta = sqrt(log(S.Overshoot/100)^2/(pi^2+log(S.Overshoot/100)^2));
Q2.zeta = 0.193;

Q2.wn =4/(Q2.zeta*S.SettlingTime);

%-------------------------

poles_G = pole(Q1.G*0.9);
impulse(Q1.G)
hold on 

pole_min = real(poles_G(2));
imp_pole = 10*pole_min;


impulse(Q1.G*41/(s+41))

Q3.p = imp_pole;
Q3.p = -41;




%-------------------------
rr = ((4+10) + (8+10))/40;
ll = (8+10) + (6+10);
cc = (1+10) + (7+10);

Q4.A = [-1/(rr*cc) 1/cc
        -1/ll -2*rr/ll];
Q4.B = [1/(rr*cc) 1/cc
        0 -rr/ll];

Q5.C = [1 rr
        -1/rr 0];

Q5.D = [0 rr  
        1/rr 0];
%run m2Submit.p;