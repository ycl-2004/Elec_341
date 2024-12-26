%%%%%%%%%%%%%%%%%%%%
% Assignment 4 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;
s = tf('s');

m1DSPlot(SN);
hold on;

L = 10*1e-3*s;
C = 1/(2*1e-6*s);
r = 50;

omega_n = 36; % 36 good 
K_dc = 1;
sigma = -9;

t = linspace(0, 0.6, 1000);
t_ms = t*1000;


ii = 1/(r+L+C);

vc = ii*C/s;

% Plot the step response

%[y, t_t] = step(vc, t_ms); % Step response with original time in seconds
%plot(t_ms, y);

%plot(t_ms, G);

s = tf('s');
G = K_dc * omega_n^2 / (s^2 - 2*sigma*s + omega_n^2);

% Plot the step response

[y, t] = step(G, t);       % Step response with original time in seconds
plot(t_ms, y);       

Q1.L = -r/(2*sigma*omega_n);
%Q1.L = 2.88888;   
%Q1.L = -r/(2*sigma*omega_n);

Q1.C = 1/(omega_n^2*Q1.L);

Rw = (4+10)/5;
Lw = ((8+10) + (8+10))/10;
Jg = ((6+10) + (1+10))/6;
Bg = ((7+10) + (8+10))/5;
Kg = (5+10)/2;

Q2.A = [0 1/Q1.C 0
        -1/(Lw+Q1.L) (-Rw-50)/(Lw+Q1.L) Kg/(Lw+Q1.L)
        0  -Kg/Jg -Bg/Jg];

Q2.B = [0
        0
        1/Jg];

C = [1 0 0
    0 1 0];

D = [0];

Q3.Gv = tf(ss(Q2.A,Q2.B,C,D));
Q3.Gi = tf(ss(Q2.A,Q2.B,C,D));
Q3.Gv = Q3.Gv(1);

Q3.Gi = Q3.Gi(2);

run m1Submit.p;


