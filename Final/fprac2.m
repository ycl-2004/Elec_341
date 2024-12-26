%%%%%%%%%%%%%%%%%%%%
% ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;
addpath('./functions/');
addpath('./practiceG/');

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

f2DSPlot(SN)
hold on;

fv = 1.83;
yline(fv);
tr = 0.0229; % value at yline first time
os = (2.063-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
w = (pi-atan(beta/zeta))/(tr*beta);

Q1.Ga = toTransfer(fv,zeta,w);

K47 = 10^5;
A78 = A*25;
B7 = B*15;
C8 = C*50;
D8 = D*E*5000;

G78 = s + A78;
H7 = (s+B7)/15;
H8 = s^2+C8*s+D8;

Q2.Hs = K47*G78/((1+H7)*(1+H8));
Q2.Ks = dcgain(Q2.Hs);
Q2.Ds=Q2.Hs/dcgain(Q2.Hs);

Z1 = A;
Z2 = 2*B;
Z3 = 7*C;
n12 = Z2/Z1;
n23 = Z3/Z2;

J1 = F*10^-4;
J2 = G*10^-3;
J3 = H*10^-2;

B2 = D;
B3 = E;
Bp = F*20;

Rw = A/2;
Lw = B*10*10^-3;
Km = (C+D)*10^-2;
Jm = G*10^-4;
Bm = (E+F)*10^-6;

Q3.Jj = J3/(n23^2)+J2+J1*n12^2+Jm*n12^2;
Q3.Bj = Bm*n12^2+B2+B3/n23^2+Bp/n23^2;
Q3.Kj = 0;

% use n12 as ratio cuz it is the ratio between
% motor and sensor?
Q4.Ye = 1/(Lw*s+Rw);
Q4.Ym = 1/(Q3.Jj*s+Q3.Bj);
Q4.Gp = feedback(Q4.Ye*Q4.Ym*n12*Km,n12*Km);

Q5.A = [-Rw/Lw -n12*Km/Lw
        n12*Km/Q3.Jj -Q3.Bj/Q3.Jj];

Q5.B = [1/Lw
        0];

Q5.C = [0 Bp/n23 
        0 1/n12];

Q5.D = [0
        0];


CF = C*100;
DC = 0.5;
NCf = D;
NCd = D+E;

deltat = 1/CF;

Q6.GHs = Q1.Ga*Q4.Gp*Q2.Hs;
poles = pole(Q6.GHs)

[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole;

Q6.Nf = (15*NCf-42)/65+DC; %3.5462
Q6.Nf = CF/(10*ceil(Q6.wd))-DC-0.5;

% WHY Add 0.5 is it Nh hold control signal or FDD
Q6.tau = (Q6.Nf+0.5)/CF; %0.002 10~25E;
%Q6.tau = 0.0022;

%Q6.tau = Q6.Nf/CF; %0.002 10~25E;
Q6.beta = exp(-deltat/Q6.tau);

%Q6.beta = exp(-1/Q6.Nf);

%---------------ASKKKKKKK-------------------------
tau2 = -Q6.tau * log(0.02);
Nterms = tau2/deltat+2; %15.XXX +2

t_values =1:1:17;
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(Nterms-1.7));
     total = total+Craw_values(i);
end

tau2 = -Q6.tau * log(0.02);
Nterms = tau2/deltat+2; %15.XXX +2

t_values =0:1:16;
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(Nterms-1.7));
     total = total+Craw_values(i);
end

Q7.W = Craw_values/total;

Q7.NC =4*Q6.tau*CF-1;
t = 0:Q7.NC-1;
Craw = @(t) exp((-4*t)/(Q7.NC-1));
result = Craw(t);
Q7.W  = result/sum(result);

Q7.NC =4*Q6.tau*CF;
Q7.num = ceil(Q7.NC);
%----------------------ASKED------------------------------

Q8.N = Q6.Nf+0.5+0.5; % 0.5 from filter 0.5 from Nh

Q8.Hc = CF/(Q8.N*s+CF)/dcgain(Q2.Hs);

Q9.G = Q1.Ga*Q4.Gp;
Q9.H = Q8.Hc*Q2.Hs;
Q9.GH = Q9.G*Q9.H;

Z1 = A;
Z2 = 2*B;
Z3 = 7*C;
%n12 = Z2/Z1;
%n23 = Z3/Z2;

loop = feedback(Q9.G,Q9.H)
Q9.Ktj = Z3/Z2;
Q9.Kjt = Z2/Z3;

%-----------------------------
s = tf('s');

Dp = CF/(s*Q8.N+CF)/s;

Q10.Dp = Dp;

[Q10.K0,pm,wx,wcp] = margin(Q10.Dp*Q9.G*Q9.H);
[K0,pm,Q10.wxo,wcp] = margin(Q10.K0*Q10.Dp*Q9.G*Q9.H);


[zz1,zz2,pm] = optimizeDoubleZero(-Q10.wxo,Q10.K0*Q10.Dp*Q9.G*Q9.H,1, deg2rad(0.5));
% JUST GUESS???????????????? -0.5+3i part JUST GUESS?????????????????
%[Zret, PMret] = newtonsCCzeroPID(-0.5 + 3i, Q10.K0 * Q10.Dp * Q9.G * Q9.H);

Q11.Z = [zz1,zz2];
Q11.PM = pm;

Dz = (s-zz1)*(s-zz2)/(zz1*zz2);
Q11.D = Dz*Q10.Dp;

step(Q11.D*Q9.G*Q9.H)
K_val = 50000000*0.875; %40.143
K_val = 50000000*0.88;  %40.0084
K_val = 100;

TotalTF = Q11.D*Q9.G*Q9.H
TargPM = 40;
K0 = 100000; % 初始值为 K = 50

% 目标函数定义
objective = @(K) abs(findPM(K, 1, TotalTF, TargPM)); % 最小化 |PM - target_PM|

% 使用 fmincon 优化
options = optimoptions('fmincon', 'Display', 'iter'); % 显示迭代过程
INT = fmincon(objective, K0, [], [], [], [], 0, [], [], options); % 限制 K > 0
Q12.K = INT; % 将计算结果赋值到 Q18.K

[~,PM,~,~] = margin(INT*Q11.D*Q9.G*Q9.H)

%K_val = 100;
%[K0,pms,wxo,wcp] = margin(K_val*Q11.D*Q9.G*Q9.H);
%Q12.K = K_val;


kpp = -1*CF/(Q8.N);
Q13.Kp  =Q12.K*(1/kpp-(zz1+zz2)/(zz1*zz2));
Q13.Ki = Q12.K;

kdp = -1*CF/(Q8.N);
p = kdp;
Q13.Kd = Q12.K*(1/(kdp)^2 - (zz1+zz2-kdp)/(kdp*(zz1*zz2)));

PID = (Q13.Kp + Q13.Ki*1/s + Q13.Kd*(-p*s)/(s-p));

Gcl = feedback(1*Q11.D*Q9.G, Q9.H);
step(Gcl)
hold on
yline(0.98);
S = stepinfo(Gcl)
Q14.Tr = S.RiseTime;
Q14.Tr = inf;
Q14.Tp = S.PeakTime;
Q14.Tp = inf;
Q14.Ts = S.SettlingTime;
Q14.Ts = inf;
Q14.OSu = S.Overshoot;
Q14.OSu = inf;
Q14.OSy = S.Overshoot;
Q14.OSy = inf;
Q14.Ess = 0;



%run f2Submit.p






