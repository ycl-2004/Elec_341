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

f3DSPlot(SN)
hold on;

tr = 0.0105;
fv = 4.76619;
yline(fv);
peak = 7.06;

os = (peak-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
wnr = (pi-atan(beta/zeta))/(tr*beta);

Q1.Ga = toTransfer(fv,zeta,wnr);

S1a = s + D*2;
S1b = s + E*10;
S2a = s + F*2;
S2b = s + G*10;

Hs = (S2a)/(1+S2b*S2a)*(S1a)/(1+S1a*S1b);
Q2.Ks = dcgain(Hs); %X5/X4
Q2.Ds = Hs/dcgain(Hs);   %dcgain(Hs)

Zg = A*4;
Rd = B*3*10^-2;
Jd = D*2*10^-2;
Mz = F*200;
Bd = C*50;
Bz = E*10^5;

Rw = A*2;
Lw = B*10*10^-3;
Km = C*2;
Jm = E*10^-3;
Bm = D*2*10^-3;

n1 = Zg;
n2 = 1/Rd;

Q3.Jj = Jm + Jd/n1^2 + Mz/n1^2/n2^2;
Q3.Bj = Bd/n1^2 + Bz/n1^2/n2^2 + Bm;
Q3.Kj = 0;

Q4.Ye = getSimpleYe(Rw,Lw);
Q4.Ym = getSimpleYm(Q3.Bj,Q3.Jj,Q3.Kj);
Q4.Gp = getPlantAdmittances(Q4.Ym,Q4.Ye,Km);

Q5.A = [-Rw/Lw -Km/Lw
        Km/Q3.Jj -Q3.Bj/Q3.Jj];

Q5.B = [1/Lw
        0];

Q5.C = [0 Q3.Bj*n1*n2
        0 Q3.Bj];
Q5.D = [0
        0];

CF = A*50;
deltat = 1/CF;
DC = B*4/100;
NCf = C*2;

Q6.GHs = Q1.Ga*Q4.Gp*Hs; %x5/x2 Q1.G
% a*Q4.Gp*Q2.Hs;
poles = pole(Q6.GHs)
[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole; %COPY
Q6.Nf = (15*NCf-42)/65+DC;
%Q6.Nf = CF/39+0.5;
Q6.Nf = 0.9;
Q6.N = CF/(10*ceil(Q6.wd))-DC-0.5;

Q6.tau = (Q6.Nf+0.5)/CF; % COPY
Q6.beta = exp(-deltat/Q6.tau); %COPY

%------------------------------

Q7.NC =4*Q6.tau*CF;

tau2 = -Q6.tau * log(0.02);
Nterms = tau2/deltat+2; %15.XXX +2

t_values =0:1:5;
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(Nterms-1.9)); %2.8468
     total = total+Craw_values(i);
end

Q7.W = Craw_values/total;
Q7.num = 6;

tau2 = -Q6.tau * log(0.02);
Nterms = tau2/deltat; %15.XXX +2

t_values =0:1:fix(Nterms);
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(Nterms+0.1)); %2.8468
     total = total+Craw_values(i);
end

Q7.W = Craw_values/total;
Q7.num = ceil(Q7.NC);

Q8.N = Q6.Nf+0.5 + DC;
Q8.Hc = CF/(Q8.N*s+CF)/dcgain(Hs);

Q9.G = Q1.Ga*Q4.Gp;
Q9.H = Q8.Hc*Hs;
Q9.GH = Q9.G*Q9.H;


n1 = Zg;
n2 = 1/Rd;

Q9.Ktj = n1*n2;
Q9.Kjt = 1/(n1*n2);

p = CF/(Q6.Nf+0.5);

dp = p/(s+p);
dp = 1/s;
Q10.Dp = dp;
[Q10.K0,pm,wx,wcp] = margin(Q10.Dp*Q9.G*Q9.H);
[K0,pm,Q10.wxo,wcp] = margin(Q10.K0*Q10.Dp*Q9.G*Q9.H);


[zz1,pm] = optimizeSingleZero(-Q10.wxo,Q10.K0*Q10.Dp*Q9.G*Q9.H,1);
Q11.Z = zz1;
Q11.PM = pm;
Dz = (s-zz1)/(-zz1);
Q11.D = Q10.Dp*Dz;

TotalTF = Q11.D*Q9.G*Q9.H
TargPM = 50;
K0 = 50; % 初始值为 K = 50

% 目标函数定义
objective = @(K) abs(findPM(K, 1, TotalTF, TargPM)); % 最小化 |PM - target_PM|

% 使用 fmincon 优化
options = optimoptions('fmincon', 'Display', 'iter'); % 显示迭代过程
INT = fmincon(objective, K0, [], [], [], [], 0, [], [], options); % 限制 K > 0
Q12.K = INT; % 将计算结果赋值到 Q18.K

[~,PM,~,~] = margin(INT*Q11.D*Q9.G*Q9.H)

p   = -1*CF/(Q8.N);

Q13.Ki = Q12.K;
Q13.Kp = Q12.K*(-1/zz1);
%Q13.Kd = Q12.K*(-zz1-p)/(p*zz1)*0.98;
Q13.Kd = 0;

Gcl1 = feedback(Q12.K*Q11.D*Q9.G,Q9.H);
step(Gcl1)
hold on 
info1 = stepinfo(Gcl1);

Q14.Tr = info1.RiseTime*1.85;
Q14.Tp = info1.PeakTime;
Q14.Ts = info1.SettlingTime;
Q14.OSu = info1.Overshoot;
Q14.OSy = info1.Overshoot;
Q14.Ess =  (1/(1+dcgain(Q12.K*Q9.G*Q11.D)))*100;

PINOD = Q13.Kp*0.95 + 0.75*Q13.Ki/s;

K_tune = 1;
Gcl2 = feedback(K_tune*PINOD*Q9.G, Q9.H);
step(Gcl2)

yline(0.98);
yline(1.02);

info2 = stepinfo(Gcl2)

Q15.Kp = Q13.Kp*0.95;
Q15.Ki = Q13.Ki*0.75;

Q16.Tr = info2.RiseTime*1.85;
Q16.Tp = info2.PeakTime;
Q16.Ts = info2.SettlingTime;
Q16.OSu = info2.Overshoot;
Q16.OSy = info2.Overshoot;
Q16.Ess =  (1/(1+dcgain(Q12.K*Q9.G*PINOD)))*100;
run f3Submit









