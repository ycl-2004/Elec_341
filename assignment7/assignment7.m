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

CF = A*10
DC = 0.3
Tf = (B+C)*10^-3
Ti = D*10^-3

CF = A*10;
DC = 30/100;
Tf = (B + C)*1e-3;
Ti = D*1e-3;
%Q1
s = tf('s')
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

Q1.G = Gp_mod;
Q1.H = Hs*Kh*Hs_gain_mod;

%Q2
Q2.N = DC+0.5
Q2.wp = CF/Q2.N
Q2.Dh = Q2.wp/(s+Q2.wp)

%Q3

Q3.NC = 4*Tf*CF+1
Q3.NC = 21;

t = 0:Q3.NC-1;
Craw = @(t) exp((-4*t)/(Q3.NC-1))
result = Craw(t)

Q3.W  = result/sum(result)
Q3.Nf = (15*Q3.NC-42)/65 *0.96
Q3.N = Q2.N + Q3.Nf
Q3.wp = CF/Q3.N



%Q4
clear s
s = tf('s')
t = 0:0.001:10;
Q4Dh = Q3.wp/(s+Q3.wp)
gain = 8.5
test = gain*Q1.G/(1+gain*Q1.G*Q1.H * Q4Dh);
[y, ~] = step(test,t)
steady_state_value = y(end)
info = stepinfo(test)
Q4.Ess = 1/(1+dcgain(gain*Q1.G*Q1.H * Q4Dh))*100
dv = (steady_state_value)+Q4.Ess/100; 
OSu = (info.Peak-dv)/dv*100
Q4.K = gain
Q4.Ts = info.SettlingTime;

%Q5
poles = pole(test)


Q5.Nf = (CF*Ti-0.5)*1.02
Q5.Nb = CF*Ti

Q5.D = 1/(Ti*15)*(s+15)/(s+1/Ti)

%Q6
Q6.Kd = 6.5*(-15+1/Ti)/(15/Ti)
gain = 6.5
test = gain*Q1.G*Q5.D/(1+gain*Q5.D*Q1.G*Q1.H * Q4Dh);
[y, ~] = step(test,t)
steady_state_value = y(end)
info = stepinfo(test)
Q6.Ess = 1/(1+dcgain(gain*Q1.G*Q1.H * Q4Dh*Q5.D))*100
dv = (steady_state_value)+Q6.Ess/100; 
OSu = (info.Peak-dv)/dv*100

Q6.Kp = gain%-Q6.Kd

Q6.Ts = info.SettlingTime

OSu

a7Submit