%%%%%%%%%%%%%%%%%%%%
% ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
addpath('./functions/');

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

SN = 48861785;
s = tf('s');

%------------------------------------------------------------
%-------------------------- Q1 ------------------------------

% Rise time copy prac2/3 
% Settle time copy prac1
x2DSPlot(SN);

hold on;

fv = 18.3;
yline(fv);
tr = 1.278*10^-3; % value at yline first time
os = (20.04-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
w = (pi-atan(beta/zeta))/(tr*beta);

Q1.Ga = toTransfer(fv,zeta,w);

%------------------------------------------------------------
%-------------------------- Q2 ------------------------------
g1 = G*100;
g2 = H*200;
G1 = g1/(s+g1);
G2 = g2/(s+g2);
Hs = -G1*G2/(G2*(-G1*G2^2-G1)-1);
Hs = Hs/s;
Hs = G1*G2/(1+G2^2*G1+(G2*G1)^2);
Q2.Ks = dcgain(Hs);
Hs = G1*G2/(1+(G2*G1)^2+(G1*G2));
Q2.Ds = Hs/dcgain(Hs)/s;

%------------------------------------------------------------
%-------------------------- Q3 ------------------------------
Dhs = A*10^-2;
Dls = B*4*10^-2;
Npw = C*10^-1; %*2*pi if needed

Rw = A*10^-2;
Lw = B*5*10^-2;
Km = C*2*10^-2;
Jm = D*10^-2;
Bm = E*10^-4;

Jhs = D*3*10^-3;
Jls = E*5*10^-2;
Jp = F*5*10^-1;
Mb = G*40;

Bhs = H*3*10^-3;
Bls = A*3*10^-3;
Bbb = B*10^-3;
Bpw = C*10^-1;
Bbw = D*20;

n1 = Dls/Dhs;
n2 = Npw*2*pi;

Q3.Mj = Jm*n1^2*n2^2 + Jhs*n1^2*n2^2 + Jls*n2^2 + Mb + Jp*n2^2;
%Q3.Mj = 2530;
%Q3.Mj = 2630;
%Q3.Mj = 2580;
%Q3.Bj = Bbb*n1^2 + Bpw*n2^2*n1^2+ Bbw*n1^2+ Bm*n1^2 + Bhs*n1^2 + Bls*n2^2*n1^2;
Q3.Bj = (Bbw*n1^2 + Bbb + Bpw*n2^2 + Bhs*n1^2 + Bls + Bm*n1^2)*1;
Q3.Bj = Bbw + (Bm + Bhs + Bbb)*n1^2*n2^2 + n2^2*(Bls+Bpw);
%Q3.Bj = 725; %728
%Q3.Bj = 714;
Q3.Kj = 0;


%------------------------------------------------------------
%-------------------------- Q4 ------------------------------

Q4.Ye = getSimpleYe(Rw,Lw);
Q4.Ym = getSimpleYm(Q3.Bj,Q3.Mj,Q3.Kj);
Q4.Gp = getPlantAdmittances(Q4.Ym,Q4.Ye,n1*n2*Km);


%------------------------------------------------------------
%-------------------------- Q5 ------------------------------

Q5.A = [-Rw/Lw -Km/Lw
        Km/Q3.Mj -Q3.Bj/Q3.Mj];

Q5.B = [1/Lw
        0];

% C matrix often deal with B and n values 
% F = B*V = B*n*w;
% torq = B*w;
Q5.C = [0 Q3.Bj*n1*n2
        0 Q3.Bj];
Q5.D = [0
        0];

CF = F*50;
DC = 0.7;
deltat = 1/CF;

Q6.GHs = Q1.Ga*Q4.Gp*Q2.Ds*Q2.Ks; %x5/x2 Q1.G
poles = pole(Q6.GHs)
[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole; %COPY
format long
poles 
Q6.wd = 0.22;
% WORK ONCE
%Q6.Nf = CF/(1800)-0.5 -DC;
Q6.Nf = CF/1800-0.5-DC; % Last version

Q6.tau = (Q6.Nf+0.5)/CF; % COPY
Q6.beta = exp(-deltat/Q6.tau); %COPY


%run x2Submit.p















