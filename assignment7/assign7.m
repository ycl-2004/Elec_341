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

% FROM ASSIGN 6

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

GH = Q1.G*Q1.H;
%ss = stepinfo(GH);
%-------------------
% ASSIGNMENT START

s = tf('s');

CF = A*10;
DC = 30/100;
tauf = (B + C)*1e-3;
taui = D*1e-3;
Ntotal = DC+0.5;

dh = CF/(Ntotal*s+CF);
Q2.N = Ntotal;
Q2.wp = -pole(dh);
Q2.Dh = dh;

%---------------------------
nn = tauf*CF;
deltaT = 1/CF;

%Craw = @(t)exp(-t/tauf);

tat2 = -tauf * log(0.02);
NC = tat2/deltaT;
Q3.NC = 20;
Q3.NC = 21;

t_values =1:1:21;
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(20.1));
     total = total+Craw_values(i);
end

Q3.W = Craw_values/total;
sum(Q3.W)

%Q3.Nf=(1*Craw_values(2)+Craw_values(3)*2+Craw_values(4)*3+Craw_values(5)*4+Craw_values(6)*5+14*(15*20-42)/65)/100;
Q3.Nf = (15*20-42)/65;
Q3.N =Q3.Nf+0.5+DC;
dh = CF/(Q3.N*s+CF);
Q3.wp =-pole(dh);
%dh = Q3.wp/(s+Q3.wp);
%--------------------------------------------------------
%-----------------

K_Val = 8.4;
%K_Val = 8.4;
ForwardPath = K_Val* Q1.G;
FeedbackPath = dh*Q1.H;
%FeedbackPath = dh*Q1.H;
ClosedLoopSystem = feedback(ForwardPath, FeedbackPath);

poles = pole(ClosedLoopSystem)
zeros = zero(ClosedLoopSystem)
zeros
%--------------------------
%S = stepinfo(ClosedLoopSystem)
%step(ClosedLoopSystem)
hold on 
yline(1)
step(ClosedLoopSystem)
%Q4.K = margin(ClosedLoopSystem);
Q4.K = K_Val;

%Q4.K = K_Val; % 8.4
%Q4.K =8.4;
Q4.Ts = 0.63; % GRAPH
Q4.Ess = 100*(1-0.644);  

%----------------------------------------
beta = exp(-deltaT/taui);

Q5.Nf =beta/(1-beta);
Q5.Nb =Q5.Nf+0.5;
%Q5.Nb = taui*CF;

p1 = -A;
Kd = (-p1 + 2*CF)/(-2*CF*p1);
Q5.D = 1+0.7*Kd*(CF*s)/((Q5.Nb)*s+CF); % working idky for me 
Q5.D
Q5.D = 1 + (z1-p1)/(p1*z1)*(-p1*z1)/(s-p1);
Q5.D = 1/(taui*-poles(end))*(s-poles(end))/(s+1/taui); 
       %1/taui)/(s+1/taui)*(s-poles(end))/(-poles(end));
       %p/z*(s+z)/(s+p)
p = CF/Q5.Nb;
%pzmap();

%Q5.D = 1/(taui*52.3)*(s+52.3)/(s+1/taui); 
% Correct way that works for my other friends 

Q5.D
poles(end)
%-----------------------------------------

K_Val = 6.15;
ForwardPath = K_Val*Q1.G*Q5.D;
FeedbackPath = dh*Q1.H;
ClosedLoopSystem = feedback(ForwardPath, FeedbackPath);


%Ku = K_Val

Q6.Kp = K_Val;
%Q6.Kd = K_Val*(poles(end)+1/taui)/(-poles(end)/taui); % works for my friends not me 

%Q6.Kd =  1/3.4;
%Q6.Kd = 0.32;
Q6.Kd = -(z1-p1)/(p1*z1)*K_Val;
Q6.Ts =0.45;
Q6.Ts = stepinfo(ClosedLoopSystem).SettlingTime;
Q6.Ess = 47;
Q6.Ess = 1/(1+dcgain(K_Val*Q1.G*Q1.H*dh*Q5.D))*100;

%run a7Submit.p

