%%%%%%%%%%%%%%%%%%%%
% Assignment 5 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;

%--------------------------------------
% Question 1 

M0 = (4+10)/5;
M1 = (8+10)/10;
M2 = (8+10)/10;
M3 = (6+10)/5;

B20 = (1+10)/2;
B21 = (7+10)/3;
B31 = (8+10)/4;

K0 = (4+10);
K1 = (8+10);
K20 = (8+10);
K32 = (6+10)/3; 



Q1.A = [-B20/M0, 0, B20/M0, 0, -1/M0, 0, 1/M0, 0 
     0, (-B21-B31)/M1, B21/M1, B31/M1, 0, -1/M1, 0, 0
     B20/M2, B21/M2, (-B20-B21)/M2, 0,0,0,-1/M2, 1/M2
     0, B31/M3, 0, -B31/M3, 0, 0, 0, -1/M3
     K0, 0, 0, 0, 0, 0, 0, 0
     0, K1, 0, 0, 0, 0, 0, 0
     -K20, 0, K20, 0, 0, 0, 0, 0
     0, 0, -K32, K32, 0, 0, 0, 0];

Q1.B = [1/M0;0;0;0;0;0;0;0];

C = [0 B31/M3 0 -B31/M3 0 0 0 -1/M3
    0 -B21 B21 0 0 0 0 0];

s = tf('s');

C = [0 0 0 1/s 0 0 0 0
    0 B21 -B21 0 0 0 0 0];
D = [0
     0];


Gd3 = ((C/(s*eye(size(Q1.A))-Q1.A))*Q1.B+D);
%Gd3 = tf(ss(Q1.A,Q1.B,C,D));

Q2.Gd3 =Gd3(1);
Q2.Gf21 =-Gd3(2);
%---------------------------------- Question 3 Begin 


Rw = (4+10)/3;
Lw = (8+10)*1e-3;
Jr = (8+10)/10*1e-3;
Br = ((6+10) + (1+10))*1e-3;

Ke = (8+10)*50*1e-3;
Km = (8+10)*50*1e-3;

Js = (4+10)*20*1e-6;
n = 5*2*pi;
Mn = (8+10)/2;
Bn = (8+10)/2;

g = 9.81;

Mtotal = M0+Mn+(Js+Jr)*n^2;

%(B20+Bn)/Mtotal 0 -B20/Mtotal 0 1/Mtotal 0 -1/Mtotal 0 Km/(Mtotal*n)

Q3.A = [(-B20-Bn-Br*n^2)/Mtotal, 0, B20/Mtotal, 0, -1/Mtotal, 0, 1/Mtotal, 0, Km/Mtotal*n
        0, (-B21-B31)/M1, B21/M1, B31/M1, 0, -1/M1, 0, 0, 0
        B20/M2, B21/M2, (-B20-B21)/M2, 0,0,0,-1/M2, 1/M2, 0 
        0, B31/M3, 0, -B31/M3, 0, 0, 0, -1/M3, 0
        K0 0 0 0 0 0 0 0 0
        0 K1 0 0 0 0 0 0 0
        -K20 0 K20 0 0 0 0 0 0
        0 0 -K32 K32 0 0 0 0 0
        -Km*n/Lw, 0, 0, 0, 0, 0, 0, 0, -Rw/Lw];

% first line -g/-10 10~ 20 error; -5> 20e -9.5 4~10e
% -8.2G ~ -8.5
Q3.B = [0 -(M0+Mn)/Mtotal*g
        0 -g
        0 -g
        0 -g
        0 0
        0 0
        0 0
        0 0
        1/Lw 0];


%---------------------------------- Question 4 Begin 

C3 = [0 0 0 0 1 0 0 0 0 
     0 0 0 0 0 1 0 0 0 
     0 0 0 0 0 0 1 0 0 
     0 0 0 0 0 0 0 1 0];

D3 = [0 0
     0 0 
     0 0
     0 0];

tfG = ((C3/(s*eye(size(Q3.A))-Q3.A))*Q3.B+D3);

Q4.Gf = tfG(3,1);
Q4.Gg = tfG(3,2);

%---------------------------------- Question 5 Begin 

supertf = tf(ss(Q3.A,Q3.B,C3,D3))*120;
supertf = (Q4.Gf)*120+Q4.Gg;

S = stepinfo(supertf);

step(supertf)
%supertf


S
Q5.Tr = S.RiseTime;
Q5.Tr = 1.22; %1.2 2~4


Q5.Tp = S.PeakTime;
Q5.Ts = S.SettlingTime;

Q5.yf = S.SettlingMax;
Q5.yf = -50;

Q5.yp = S.Peak;
Q5.yp = S.SettlingMin; 

Q5.OSy = S.Overshoot;

run a5Submit.p;
