%-------------- Question 1-------------------
%--------------------------------------------
% In underdamp-condition, Check which time type they asked
% 
% Settle Time, 當graph的值到達FV的98%
fv = 4.764;
st = 0.086;
peak = 6;

os = (peak-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
wns = 4/(st*zeta);
Q1.Ga = toTransfer(fv,zeta,wns);

function [f] = toTransfer(FV, ZETA, Wn)
    num = FV*Wn^2;
    den = [1 2*ZETA*Wn Wn^2];
    f = tf(num,den);
end

% Rise Time, 當graph的值第一次到達FV

tr = 0.0229;
fv = 1.83;
peak = 2.063;

os = (peak-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
wnr = (pi-atan(beta/zeta))/(tr*beta);

Q1.Ga = toTransfer(fv,zeta,wnr);

% Peak Time, 當graph的值第一次到達Peak Value 

tp = 1;
fv = 1;
peak = 1;

os = (peak-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
wnp = pi/(tp*beta);

Q1.Ga = toTransfer(fv,zeta,wnp);

% In Overdamp-condition, ----------------------------------
% Find Settle Time, time it reaches FV (Works for assign2)
% OR taut, time it takes to 63% of FV

% tr1 
st = 1;
taut = st/4;

Tat90FV = 2;
Tat10FV = 1;

tr1 = Tat90FV - Tat10FV;
wn5 = pi/0.63; % TEST values
zeta = 1/(wn5*taut);


%-------------- Question 2-------------------
%---------------TF Function------------------
%----------USUALLY SHOULD BE FIND------------
%------------CHAT GPT IF NEEDED--------------



%-------------- Question 3-------------------
%--------Pinion -> Smaller of Gear ----------
%-----------Gear -> Larger of Gear ----------
%------------N = GearN/PinionN --------------
% 兩個齒輪 找 N ratio by 2pi*r1/(2pi*r2), N  = Z2/Z1 = largeZ/smallZ

% r1 -> small one then = 1
%Its from the arc length formula
%D=rtheta
%Where D is the arc length
%We want to find theta/D for the gear ratio of a pulley
%So it’s 1/r

Q3.Jj = Jm + Jd/n1^2 + Mz/n1^2/n2^2;
Q3.Bj = Bd/n1^2 + Bz/n1^2/n2^2 + Bm;
Q3.Kj = 0;


%-------------- Question 4-------------------
%-------------Ideally All copy---------------

Q4.Ye = getSimpleYe(Rw,Lw);
Q4.Ym = getSimpleYm(Q3.Bj,Q3.Jj,Q3.Kj);
Q4.Gp = getPlantAdmittances(Q4.Ym,Q4.Ye,Km)/s; % check for 1/s is it in Gp

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

%-------------- Question 5-------------------
%-------------Ideally All copy---------------
%--------------------------------------------
% Adjust the -Km terms with n if applied *n^2 in part 3

% A and B have [iw
%               wj] and [x3] as the state vector 

Q5.A = [-Rw/Lw -Km/Lw
        Km/Q3.Jj -Q3.Bj/Q3.Jj];

Q5.B = [1/Lw
        0];

Q5.C = [1 1
        1 1];
Q5.D = [0
        0];


%-------------- Question 6-------------------
%-------------Ideally All copy---------------
%--------------------------------------------
% Adjust Nf mostly only 

Q6.GHs = Q1.Ga*Q4.Gp*Q2.Hs; %x5/x2 Q1.G
poles = pole(Q6.GHs)
[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole; %COPY

Q6.Nf = (15*NCf-42)/65+DC; % WORK ONCE
Q6.Nf = 0.9;
Q6.Nf = CF/(10*ceil(Q6.wd))-0.5;
Q6.Nf = CF/(Q6.wd*10) -DC - 0.5; % Last version

Q6.tau = (Q6.Nf+0.5)/CF; % COPY
Q6.beta = exp(-deltat/Q6.tau); %COPY

%-------------- Question 7-------------------
%-------------Ideally All copy---------------
%--------------------------------------------
% Adjust -1.9 terms by 0.1 and 

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


Q7.NC =4*Q6.tau*CF+1;
t = 0:Q7.NC-1;
Craw = @(t) exp((-4*t)/(Q7.NC-1));
result = Craw(t);
Q7.W  = result/sum(result);


Q7.NC =4*Q6.tau*CF;
Q7.num = ceil(Q7.NC);


%-------------- Question 8-------------------
%-------------Ideally All copy---------------
%--------------------------------------------
% HOPE NO ADJUSTMENT NEEDED SHOULD PASS HERE

Q8.N = Q6.Nf+0.5 + DC;
Q8.Hc = CF/(Q8.N*s+CF)/dcgain(Q2.Hs);

%-------------- Question 9-------------------
%-------------Ideally All copy---------------
%--------------------------------------------
% HOPE NO ADJUSTMENT NEEDED 

Q9.G = Q1.Ga*Q4.Gp;
Q9.H = Q8.Hc*Q2.Hs;
Q9.GH = Q9.G*Q9.H;

Q9.Ktj = 1;
Q9.Kjt = 1; % Play with the n ratio and see;
            % dont cry dont cry dont cry 

 
%-----------------PID------------------------
%-------------- Question 10------------------
%-------------Ideally All copy---------------
%--------------------------------------------

% -> ON NOTES, 2*CF/(s*())
Dp = CF/(s*Q8.N+CF)/s;
Q10.Dp = Dp;

%------------------PI-----------------------
Dp = 1/s;
%-------------------------------------------


%------------------PD-----------------------
p_Dp = -CF/Nt; %Nt=0.5 prob pr Q6.Nf + 0.5;
Dp = -p_Dp/(s-p_Dp);
%-------------------------------------------

[Q10.K0,pm,wx,wcp] = margin(Q10.Dp*Q9.G*Q9.H);
[K0,pm,Q10.wxo,wcp] = margin(Q10.K0*Q10.Dp*Q9.G*Q9.H);

%-------------- Question 11------------------
%-------------Ideally All copy---------------
%-------------Zero PM Control----------------

% -> If requirement, (WnRes = 1, ZetaRes = 0.01) -> (1,deg2rad(0.5))
%                                                ->  Work for F2
%                    (WnRes = 0.1, ZetaRes = 0.01) -> (0.1,deg2rad(1))
%                                                  -> Work for Pro2

[zz1,zz2,pm] = optimizeDoubleZero(-Q10.wxo,Q10.K0*Q10.Dp*Q9.G*Q9.H,1, deg2rad(0.5));
Q11.Z = [zz1,zz2];
Q11.PM = pm;
Dz = (s-zz1)*(s-zz2)/(zz1*zz2);
Q11.D = Dz*Q10.Dp;

[zz1,pm] = optimizeSingleZero(-Q10.wxo,Q10.K0*Q10.Dp*Q9.G*Q9.H,1);
Q11.Z = zz1;
Q11.PM = pm;
Dz = (s-zz1)/(-zz1);
Q11.D = Q10.Dp*Dz;


%-------------- Question 12------------------
%-------------Ideally All copy---------------
%-------------Adjust K to get pm-------------
K_val = 50000000*0.875;
[K0,pms,wxo,wcp] = margin(K_val*Q11.D*Q9.G*Q9.H);
Q12.K = K_val;


TotalTF = Q11.D*Q9.G*Q9.H
TargPM = 40;
K0 = 50; % 初始值为 K = 50

% 目标函数定义
objective = @(K) abs(findPM(K, 1, TotalTF, TargPM)); % 最小化 |PM - target_PM|

% 使用 fmincon 优化
options = optimoptions('fmincon', 'Display', 'iter'); % 显示迭代过程
INT = fmincon(objective, K0, [], [], [], [], 0, [], [], options); % 限制 K > 0
Q12.K = INT; % 将计算结果赋值到 Q18.K

[~,PM,~,~] = margin(INT*Q11.D*Q9.G*Q9.H)


%-------------- Question 13------------------
%-------------Ideally All copy---------------
%-------------Adjust p to get Ks-------------
p   = -1*CF/(Q8.N);
kpp = -1*CF/(Q8.N);
kdp = -1*CF/(Q8.N);

Q13.Kp  =Q12.K*(1/kpp-(zz1+zz2)/(zz1*zz2));
Q13.Ki = Q12.K;
Q13.Kd = Q12.K*(1/(kdp)^2 - (zz1+zz2-kdp)/(kdp*(zz1*zz2)));



% FINDING Ts
step(Gcl)
hold on

yline(1.02)
yline(0.98)


%------------------------------------------------------------
%-------------------------- Q6 ------------------------------

CF = 1;
DC = 1;

Q6.GHs = Q1.Ga*Q4.Gp*Hs; %x5/x2 Q1.G
poles = pole(Q6.GHs)
[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole; %COPY

Q6.Nf = (15*NCf-42)/65+DC; % WORK ONCE
Q6.Nf = CF/(10*ceil(Q6.wd))-0.5 -DC;
Q6.Nf = CF/(Q6.wd*10) -DC - 0.5; % Last version

Q6.tau = (Q6.Nf+0.5)/CF; % COPY
Q6.beta = exp(-deltat/Q6.tau); %COPY

%------------------------------------------------------------
%-------------------------- Q7 ------------------------------

Q7.NC =4*Q6.tau*CF;
Q7.num = ceil(Q7.NC);

%------------------------------------------------------------
%-------------------------- Q8 ------------------------------

Q8.N = Q6.Nf+0.5 + DC;
Q8.Hc = CF/(Q8.N*s+CF)/dcgain(Hs);

%------------------------------------------------------------
%-------------------------- Q9 ------------------------------

Q9.G = Q1.Ga*Q4.Gp;
Q9.H = Q8.Hc*Hs;
Q9.GH = Q9.G*Q9.H;

Q9.Ktj = 1;
Q9.Kjt = 1;

%------------------------------------------------------------
%-------------------------- Q10 ------------------------------