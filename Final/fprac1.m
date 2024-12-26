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

f1DSPlot(SN)
hold on 

fv = 4.764;

yline(fv*0.98);
yline(fv*1.02);

st = 0.086;
os = (6-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
wn = 4/(st*zeta);

Q1.Ga = toTransfer(fv,zeta,wn);

%-------------------------------------

Ko = 3*10^6;
p1 = B*25;
p2 = C*25;

h1 = 1/(s+p1);
h2 = 1/(s+p2);

Hs = (Ko*h1*h2)/(1+h1+h2)
Q2.Ks = dcgain(Hs);
Q2.Ds = Hs/dcgain(Hs);


%-------------------------------------

Dp = A*0.01;
%Dpr = Dp/2;
Jp = B*2*10^-4;
Mt = C/2;
Bp = D*10^-5;
Bw = E + F;
Kt = G/2;

Rw = A/3;
Lw = B/5;
Km = C*2*10^-2;
Jm = D*4*10^-4;
Bm = E*2*10^-6;

n = (2/Dp); % n = 2/Dp

Q3.Jj = Jp + Mt/n^2+Jm;
Q3.Bj = Bp + Bw/n^2+Bm;
Q3.Kj = Kt/n^2;

Q4.Ye = getSimpleYe(Rw,Lw);
Q4.Ym = getSimpleYm(Q3.Bj,Q3.Jj,Q3.Kj);
Q4.Gp = getPlantAdmittances(Q4.Ym,Q4.Ye,Km)/s;

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

Q5.A = [-Rw/Lw -Km/Lw 0 
        Km/Q3.Jj -Q3.Bj/Q3.Jj -1/Q3.Jj
        0 Q3.Kj 0];

Q5.B = [1/Lw
        0
        0];

Q5.C = [0 Bw/n 0
        0 0 n];
Q5.D = [0
        0];

CF = D*20;
DC = (E+F)/100;
deltat = 1/CF;

% NOTIES -> Q2.Ds*Q2.Ks = Hs
Q6.GHs = Q1.Ga*Q4.Gp*Q2.Ds*Q2.Ks; %x5/x2 Q1.G 
poles = pole(Q6.GHs)
[~, idx] = max(real(poles));
dominant_pole = poles(idx);
Q6.wd = -dominant_pole; %COPY
Q6.wd = 0.28;

Q6.Nf = CF/(10*ceil(Q6.wd))-DC-0.5;

Q6.tau = (Q6.Nf+0.5)/CF; % COPY
Q6.beta = exp(-deltat/Q6.tau); %COPY

Q7.NC =4*Q6.tau*CF+1;
t = 0:Q7.NC-1;
Craw = @(t) exp((-4*t)/(Q7.NC-1));
result = Craw(t);
Q7.val  = result/sum(result);
Q7.num = 128;

Q8.N = Q6.Nf+0.5 + DC;
Q8.Hc = CF/(Q8.N*s+CF)/dcgain(Hs);

Q9.G = Q1.Ga*Q4.Gp;
Q9.H = Q8.Hc*Hs;
Q9.GH = Q9.G*Q9.H;

Q9.Ktj = n;
Q9.Kjt = 1/n; % Play with the n ratio and see;
            % dont cry dont cry dont cry 

Dp = CF/(s*Q8.N+CF)/s;
Q10.Dp = Dp;
[Q10.K0,pm,wx,wcp] = margin(Q10.Dp*Q9.G*Q9.H);
[K0,pm,Q10.wxo,wcp] = margin(Q10.K0*Q10.Dp*Q9.G*Q9.H);

%WnRes = 0.1;
%ZetaRes = 0.01;

%Open_Loop2 = Q10.K0*Q10.Dp*Q9.G*Q9.H;
%wnz_range = 0:WnRes:5; % Range for wnz
%zeta_range = 0:ZetaRes:2; % Range for zeta

% Pre-allocate a matrix to store phase margin values
%PM_values = zeros(length(zeta_range), length(wnz_range));

% Loop through all combinations of wnz and zeta


%for i = 1:length(wnz_range)
%    for j = 1:length(zeta_range)
%        wnz = wnz_range(i);
%        zeta = zeta_range(j);

        % Define the transfer function for the current wnz and zeta
%        Q17_OL = Open_Loop2 * (s^2 + 2*zeta*wnz*s+ wnz^2)/wnz^2;

        % Calculate the phase margin
%        [~, PM, ~, ~] = margin(Q17_OL);

        % Store the phase margin value
%       PM_values(j, i) = PM;
%    end
%end

% Create a 3D surface plot

%[Wnz, Zeta] = meshgrid(wnz_range, zeta_range) % Create a grid for plotting
%figure
%surf(Wnz, Zeta, PM_values, 'EdgeColor', 'none') % 3D surface plot

% Customize the plot
%xlabel('\omega_{nz}', 'Interpreter', 'tex');
%ylabel('\zeta', 'Interpreter', 'tex');
%zlabel('Phase Margin (degrees)', 'Interpreter', 'tex');
%title('Phase Margin vs. \omega_{nz} and \zeta');
%colorbar; % Add a color bar for reference
%grid on;
%view(45, 30); % Adjust the viewing angle


%[Zret, PMret] = newtonsCCzeroPID(0.3,0.5, Open_Loop2);
%[Zret, PMret] = newtonsCCzeroPID(0.55,0.7, Open_Loop2);
%Q11.Z = [Zret conj(Zret)]; %transpose(roots([1 2*zetaopt*wnzopt wnzopt^2]));

Q11.Z = [-0.3217+0.3678i -0.3217-0.3678i];
%Q11.PM = PMret; %;maxPM; %94.792;
Q11.PM = 94.792;
Q11.D = Q10.Dp/(Q11.Z(1) * Q11.Z(2)) * (s-Q11.Z(1))*(s-(Q11.Z(2)));
%step(Q11.D)

K_val = -30.1;
[K0,pms,wxo,wcp] = margin(K_val*Q11.D*Q9.G*Q9.H);
Q12.K = K_val;


%TotalTF = Q11.D*Q9.G*Q9.H
%TargPM = 60;
%K0 = 50;


%objective = @(K) abs(findPM(K, 1, TotalTF, TargPM));

%options = optimoptions('fmincon', 'Display', 'iter'); 
%INT = fmincon(objective, K0, [], [], [], [], 0, [], [], options); 
%Q12.K = INT;
Q12.K = 0.1857;

%[~,PM,~,~] = margin(INT*Q11.D*Q9.G*Q9.H)

zz1 = Q11.Z(1);
zz2 = Q11.Z(2);

p   = -1*CF/(Q8.N);
kpp = -1*CF/(Q8.N);
kdp = -1*CF/(Q8.N);

Q13.Kp  =Q12.K*(1/kpp-(zz1+zz2)/(zz1*zz2));
Q13.Ki = Q12.K;
Q13.Kd = Q12.K*(1/(kdp)^2 - (zz1+zz2-kdp)/(kdp*(zz1*zz2)));

Gcl = feedback(Q12.K*Q9.G*Q11.D,Q9.H)
%step(Gcl)
%hold on 
info = stepinfo(Gcl,'RiseTimeLimits',[0 1])

Q14.Tr = info.RiseTime;
Q14.Tp = info.PeakTime;
Q14.Ts = info.SettlingTime;
Q14.OSu = info.Overshoot;
Q14.OSy = info.Overshoot;
Q14.Ess = (1/(1+dcgain(Q12.K*Q9.G*Q11.D)))*100;

Kin.K = 1;
Kin.Kd = Q13.Kd/2;
Kin.Ki = Q13.Ki/1.4;
Kin.Kp = Q13.Kp;

PID = Kin.K*(Kin.Kp + Kin.Ki*1/s + Kin.Kd*(-p*s)/(s-p));
K_tune = 0.99; %0.284 NEXT TIME TRY 0.6 up by 0.01
Gcl2 = feedback(K_tune*PID*Q9.G, Q9.H);
%step(Gcl2)
info2 = stepinfo(Gcl2,'RiseTimeLimits',[0,1]);
hold on

%heurRCGTune(Kin,Q14.Ts*0.5,Q14.OSu*0.7/100,0,p,Q9.G,Q9.H);

Q15.K = 1;
Q15.Kp = Kin.Kp;
Q15.Ki = Kin.Ki;
Q15.Kd = Kin.Kd;

Q16.Tr = info2.RiseTime*0.98;
Q16.Tp = info2.PeakTime;
Q16.Ts = info2.SettlingTime;
Q16.OSu = info2.Overshoot;
Q16.OSy = info2.Overshoot;
Q16.Ess = (1/(1+dcgain(Q15.K*Q9.G*PID)))*100;

run f1Submit.p;





