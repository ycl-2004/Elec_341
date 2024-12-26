%%%%%%%%%%%%%%%%%%%%
% Assignment 3 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;

%---------------------------
%Question 1 Parameters

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

%---------------------------
%Question 1 

R20 = 1/B20;
R21 = 1/B21;
R31 = 1/B31;

L0 = 1/K0;
L1 = 1/K1;
L20 = 1/K20;
L32 = 1/K32;

C0 = M0;
C1 = M1;
C2 = M2;
C3 = M3;

F0 = 1;

s = tf('s');

Y = [1/R20+1/(s*L20)+1/(s*L0)+s*C0 0 -1/(s*L20)-1/R20 0
    0 1/(s*L1)+1/R21+1/R31+s*C1 -1/R21 -1/R31
    -1/(s*L20)-1/R20 -1/R21 1/(s*L32)+1/R21+1/R20+1/(s*L20)+s*C2 -1/(s*L32)
    0 -1/R31 -1/(s*L32) 1/R31+1/(s*L32)+s*C3];

I = [F0; 0; 0; 0];
V = Y\I;
V = V/s;

Q1.Gd1 = V(2);
Q1.Gd3 = V(4);

Q1.Gd1
Q1.Gd3

% Question 1 End
%---------------------------
% Question 2

Q2.Gf1 = -1*V(2)*K1;
Q2.Gf32 = -1*(V(4)-V(3))*K32;

% Question 2 End
%---------------------------

%---------------------------
%Question 3 4 Parameters

Rw = (4+10)/3;
Lw = (8+10)*1e-3;
Jr = (8+10)/10*1e-3;
Br = ((6+10) + (1+10))*1e-3;

Ke = (8+10)*50*1e-3;
Km = (8+10)*50*1e-3;


%----------------------------
% Question 3 
% Vin = d_Iw*Lw + Iw*Rw + Ke*W
% Km*Iw = W*Jr*s + W*Br
Q3.Gi = 1/((s*Lw)+Rw+(Ke*Km)/(Jr*s+Br));


% Question 3 End 
%-----------------------------
% Question 4

clear s;

syms w iw s vin;

eq1 = iw == (Jr*s*w+Br*w)/(Km);
eq2 = vin == Lw*s*iw+iw*Rw+Ke*w;

[vIn, wIn] = solve(eq1, eq2, vin, w);
vIn

tf4 = wIn/vIn;
tf4 = simplify(tf4);

[num, den] = numden(tf4);

clear s;

s = tf("s");

t4f = tf(sym2poly(num), sym2poly(den));
t4f
Q4.Gw = t4f;

% Question 4 End
%----------------------------------
% Question 5
t_val = 0:1e-3:50e-3;

current = Q3.Gi/s;
speed = Q4.Gw/s;

syms s t
[symNum, symDen] = tfdata(current, 'v'); 
current_syms = poly2sym(symNum, s) / poly2sym(symDen, s);
current_t = ilaplace(current_syms, s, t);

[symNum, symDen] = tfdata(speed, 'v'); 
omega_syms = poly2sym(symNum, s) / poly2sym(symDen, s);
omega_t = ilaplace(omega_syms, s, t);

current_vals = subs(current_t, t, t_val);
omega_vals = subs(omega_t,t,t_val);

current_vals = double(current_vals);
omega_vals = double(omega_vals);

current_Rw = current_vals.^2;
p_Rw = current_Rw * Rw;
p_Br = omega_vals.^2 * Br;

Q5.Ed = sum(p_Rw + p_Br) * 1e-3;

E_Lw = 0.5 * Lw * current_vals(end)^2; 
E_Jr = 0.5 * Jr * omega_vals(end)^2;  
Q5.Es = E_Lw + E_Jr;

Q5.Ep = sum(current_vals) * 1e-3;

run a3Submit.p;
