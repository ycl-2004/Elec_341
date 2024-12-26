% Second order approx notes 
% FV -> Final Value 
% Osy -> 
% Tr -> Rise time, min time for step response reach value final value 
% Tp -> Peak time, max value of step function, time that occur
% Ts -> Settle time, step response settle to 98% of FV


%  OSy = (Peak - FV)/FV*100;
%  OSu = (Peak-1)*100;
%  Ess = 1/(1+dcgain(K*D*G*H))*100;

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

%---- FINAL USE FUNCTION -----------------------

function [f] = toTransfer(FV, ZETA, Wn)
    num = FV*Wn^2;
    den = [1 2*ZETA*Wn Wn^2];
    f = tf(num,den);
end

function Ye = getSimpleYe(Rw,Lw)
    s = tf('s');
    Ye = 1/(s*Lw + Rw);
end

function Ym = getSimpleYm(Bt,Jt,Kt)
    s = tf('s');
    Ym = 1/(Bt + s*Jt + Kt/s);
end

function Gp = getPlantAdmittances(Ym,Ye,Km, Ga)
    if ~exist('Ga','var')
        Ga = 1;
    end
    Gp = Ga*feedback(Ye*Km*Ym, Km);
end

%----------------------------------------------------------

clear all; clc;
addpath '..'
save_student_number

s = tf('s');


CF = C_sn * 30;
DC = D_sn * 4; % percent
Nw = E_sn + F_sn;

Rw = A_sn;
Lw = B_sn * 10 / 1000;
Km = C_sn + D_sn;
Bm = (E_sn+F_sn)*1e-4;
Jm = G_sn*1e-4;

S1a = s + A_sn + B_sn;
S1b = s + (C_sn + D_sn)*5;
S2a = s + (E_sn + F_sn);
S2b = s + (G_sn + H_sn)*5;

N = A_sn + B_sn + C_sn + D_sn;
Rd = (E_sn + F_sn + G_sn)*1e-2;
N2 = 1/Rd;
Bd = A_sn * 50;
Jd = (E_sn + F_sn)*1e-2;
Bz = (A_sn + B_sn)*100;
Mz = (G_sn + H_sn)*100;


x2DSPlot(SN)
[t,v] = getFigData(1);
Tp = 0.01326;
tau = getTimeConstant(t,v) / 1000;
Tr = (getRiseTime(t,v) +3.9)/ 1000 ;
plotSettleLines(v)
Ts = 24.6/1000;
peak = 5.91038;
OS = peak/4.43333-1;
[wnp, wnr, wns, zeta, beta] = getOmegaNs(Tp,Tr,Ts,OS);
wn = wnr;
s = tf('s');
Ga = wn^2/(s^2+2*zeta*wn*s + wn^2) * 4.43333;
Q1.Ga = Ga;

Hs1 = S1a/(1 + S1a*S1b);
Hs2 = S2a/(1+S2a*S2b);
Hs = Hs1*Hs2;
Ks = dcgain(Hs);

Q2.Hs1 = Hs1;
Q2.Hs2 = Hs2;
Q2.Hs = Hs;
Q2.Ks = Ks;

Bd_ = Bd/N^2;
Bz_ = Bz/N^2/N2^2;
Jd_ = Jd/N^2;
Mz_ = Mz/N^2/N2^2;

Bt = Bm+Bd_+Bz_;
Jt = Jm + Jd_ + Mz_;

Q3.B = Bt;
Q3.J = Jt;

Ye = getSimpleYe(Rw, Lw);
Ym = getSimpleYm(Bt, Jt, 0);
Gp = getPlantAdmittances(Ym,Ye,Km,1);

Q4.Ye = Ye;
Q4.Ym = Ym;
Q4.Gp = Gp;

A = [[-Rw/Lw, -Km/Lw]; [Km/Jt, -Bt/Jt]];

Q5.A = A;

Kf = 1/Q2.Ks/N/N2;
Ncalc = DC/100;

[Nws, cs] = generateWSCoefficients(Nw);

Nt = Nws + 0.5 + Ncalc;

Dh = CF/(Nt*s + CF);

Q6.Kh = Kf;
Q6.Nf = Nws;
Q6.Nt = Nt;
Q6.Dh = Dh;

Q7.C = [cs(1), cs(2), cs(3), cs(4), cs(5)];

G = Gp*Ga;
H = Hs*Kf*Dh;
GH = G*H;
Kz = 1/N/N2;
Ku = 1.73e3;
XS = 7;

Q8.G = G;
Q8.H = H;
Q8.GH = GH;
Q8.Kz = Kz;
Q8.Ku = Ku;
Q8.XS = XS;


Ntdd = 0.5 + Nws;
invNd = 1/Ntdd;
p = -1*invNd*CF;

Dp = -p/(s-p);
Gol = Dp*GH;
K0 = margin(Dp*GH);
[gm, pm, wxo] = margin(K0*Gol);

Q9.Nt = Ntdd;
Q9.Dp = Dp;
Q9.K0 = K0;
Q9.wxo = wxo;
pdy = K0 * Dp * GH;

[oz, PM] = optimizeSingleZero(-wxo, pdy, 1);
% oz = -wxo-81;
% oz = -wxo;
% oz = -3.4083;

Dz = (s-oz)/(-oz);


D = Dp*Dz;
PM = phsMargin(K0*D*GH);
Q10.Z = oz;
Q10.PM = PM;
Q10.D = D;

K = K0/2.1;
Kp = -oz;
Kd = K*(1/p - 1/oz);
Q11.K = K;
Q11.Kp =8.883e06/8052;
Q11.Kd = Kd;

% PD = K_mod*(Kp_mod  + Kd*(-p*s)/(s-p));
    
Q12.Ess = 43.6874;
%x2Submit


wnz = 0;
zeta = 0.01;
maxPM = -inf;
while wnz < 5
    while zeta < 5
        [~, PM, ~, ~] = margin(Q10.K0*Q10.Dp*Q9.H*Q9.G * (s^2 + 2*zeta*wnz*s+ wnz^2)/wnz^2);
        if maxPM < PM 
            maxPM = PM;
            zetaopt = zeta;
            wnzopt = wnz;
        end
    fprintf('%.2f %.2f\n', wnz, zeta);
    zeta = zeta + 0.01;
    end
    wnz = wnz + 1;
    zeta = 0.01;
end

wnzopt = wnzopt+1
[~, maxPM, ~, ~] = margin(Q10.K0*Q10.Dp*Q9.H*Q9.G * (s^2 + 2*zetaopt*wnzopt*s+ wnzopt^2)/wnzopt^2)
something = roots([1 2*zetaopt*wnzopt wnzopt^2]);
Q11.Z = [something(1) something(2)]; 
Q11.PM = maxPM;
Q11.D = Q10.Dp/(Q11.Z(1) * Q11.Z(2)) * (s-Q11.Z(1))*(s-(Q11.Z(2)));

zeta = -200;
maxPM = -inf;

while zeta < 200
        [~, PM, ~, ~] = margin(Q10.K0Q10.DpQ9.HQ9.G (s+zeta)/(zeta));

        if maxPM < PM 
            poles = pole(feedback(Q10.K0Q10.DpQ9.G *(s+zeta)/(zeta),Q9.H))
            if all(real(poles) < 0)
                maxPM = PM;
                zetaopt = zeta;
            end
        end
        fprintf('%.2f\n', zeta);
        zeta = zeta + MagRes;
end