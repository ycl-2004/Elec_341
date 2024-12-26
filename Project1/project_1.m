%%%%%%%%%%%%%%%%%%%%
% Project_1 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%
clear all;
clc



SN = 48861785;


SN = 48861785;
p1DSPlot(SN);
%------------------------Question 1 Amplifier Metrics
% FV -> Final Value
% OS -> Overshoot (Peak-FV)/FV *100%
% Tr -> Rise_Time, time it takes to the final value
% Tp -> Peak_Time, time it takes to reach the first maximum overshoot
% Ts -> Settle_Time, time it takes for system to stable

Q1.Tr = 0.00073; % 0.00072~74 G
Q1.Tp = 0.001705; % 0.00168~173 G
Q1.Ts = 0.004695; % 0.00460~479 G

FV1 = 18;
overshoot1 = (22.5-FV1)/FV1; % G
Q1.OSy = overshoot1*100;

%-----------------------Question 2 Underdamp
zeta2 = -log(overshoot1) / sqrt(pi^2 + (log(overshoot1))^2); 
beta2 = sqrt(1-zeta2^2);

wn2 = (pi-atan(beta2/zeta2))/(Q1.Tr*beta2);

Q2.Ga =toTransfer(FV1,zeta2,wn2);

%-----------------------Question 3 Transfer Function

syms U Y X s A;

P1 = (4+10)*7;
P2 = (8+10)*800;
P3 = (8+10)*8;
P4 = (6+10)*700;
P5 = (1+10)*600;
P6 = (7+10)*50;
P7 = (8+10)*500;
P8 = (5+10)*5;

G1 = P1/(s+P2);
G2 = P3/(s+P4);
G3 = 10^5/(s+P5);
G4 = P6/(s+P7);

H1 = 4/(s+P8);

%eq1 = (U*G1 + X)*G2 == Y;
%eq2 = (U - (Y+X*G4)*H1)*G3 == X;

%[X, Y] = solve(eq1, eq2, X, Y);

%Hs = (G3+G1*G4*H1*G3+G1)/(H1*G3*(1+G4/G2)+1/G2);
%Hs = simplify(Hs);

eq3 = (U-(Y+A*G4)*H1)*G3 == A;
aa = solve(eq3,A);
A = aa;

eq4 = (U*G1+A)*G2 == Y;
 
tftftf = solve(eq4,Y);

tftftf = tftftf/U;

[num, den] = numden(tftftf);

t3 = tf(sym2poly(num), sym2poly(den));


%[num, den] = numden(Hs);

%t3 = tf(sym2poly(num), sym2poly(den));

Q3.Ks = dcgain(t3)*0.001;

Q3.Ds = t3 /dcgain(t3);

clear s;
clear A;

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

%------------------------Question 4
s = tf("s");

Rw = (A)/2; %ohm
Lw = (B)*30*1e-6; % H
Km = (C)*1e-3; % Nm/A
Jr = (D)/15*1e-7; %g-cm^2 ->Kg-m^2
Br = (E)/30*1e-6; %Nms
Mm = ((F)+(G)); %kg

Js = (A)/5*1e-7; %g-cm^2 -> Kg-m^2
Ms = (B)/4*1e-3; %g -> Kg     
Bs = (C)/3;      %Ns/m
Ns = (D)/(2*E)*100*(2*pi); %turn/cm -> rad/m;
Jf = F/3*1e-7; %g-cm^2 -> Kg-m^2
Bf = G*1e-3;   %Nms
Nf = (H)*5*pi/9; % deg/cm -> rad/m
Bt = (A)*1e-3; %Nms
Kt = (B)*(C)*1e-3; %Nm

Bl = (D)/5; %Ns/m
Kl = (E)*(F); %N/m
L6 = ((G)+(H))*4*1e-3; %mm -> m


Nt = L6;
%Nt = Bt/Bf;

Q4.Ye = minreal(1/(s*Lw + Rw));


%------------------------Question 5

Jtotal = Jr + Js + (Mm+Ms)/Ns^2 + Jf/Ns^2*Nf^2;
Btotal = Br + Bs/Ns^2 + Bf*Nf^2/Ns^2 + Bt*Nf^2/Ns^2 + Bl*Nf^2*Nt^2/Ns^2;
Ktotal = Kt*Nf^2/Ns^2+Kl*Nt^2*Nf^2/Ns^2;

%Q5.Ym = Km/((Rw+Lw*s)*((Jr+Mm)*s+Br)+Km^2);
Q5.Ym = minreal(1/(s*Jtotal + Btotal + Ktotal/s));
%Q5.Ym = Km^2/((Rw+Lw*s)*(Jr+Mm)*s+Br*(Rw+Lw*s)+Km^2);

ratio = 3;

Q5.Ym = 1/(Jr*s+Br+Js*s+(Mm*s+Ms*s+Bs+Nf^2*(ratio*Jf*s+ratio*Bf+RR(ratio*Kt/s+ratio*Bt,Nt^2*(ratio*Kl/s+ratio*Bl))))/Ns^2);

%------------------------Question 7
Q6.taug = (Mm)*9.81/Ns;

Q7.Gtau = Km/(Rw+s*Lw);


%------------------------ Copy 

J1 = Jr + Js;

M1 = Mm + Ms;

Jt = J1 + Ns^2*M1 + Ns^2*Nf^2*Jf;

Bf_ = Bf*Ns^2*Nf^2;
Bs_ = Bs*Ns^2;

Ba = Br + Bf_ + Bs_;

Kt_ = Kt*Ns^2*Nf^2;
Bt_ = Bt*Ns^2*Nf^2;
Kl_ = Kl*Ns^2*Nf^2*Nt^2;
Bl_ = Bl*Ns^2*Nf^2*Nt^2;

%%%%%


Ra = 1/Ba;
Rt = 1/Bt_;
Rl = 1/Bl_;
Ct = Jt;
Lt = 1/Kt_;
Ll = 1/Kl_;

zRa = Ra;
zRt = Rt;
zRl = Rl;
zCt = 1/(s*Ct);
zLt = s*Lt;
zLl = s*Ll;

yRa = 1/zRa;
yRt = 1/zRt;
yRl = 1/zRl;
yCt = 1/zCt;
yLt = 1/zLt;
yLl = 1/zLl;

Y = [yRa + yCt + yRt + yLt, -yRt-yLt;
     -yRt-yLt, yRt+yLt+yRl+yLl];
I_vec = [1;0];
V = Y\I_vec;

beta1 = -Kl_/(Bl_+Bt_);
beta2 = -1/(Bl_+Bt_);
beta3 = -Bt_/(Bl_+Bt_);

b1 = (Bt_*(beta3-1)-Ba)/Jt;
b2 = (1+beta2*Bt_)/Jt;
b3 = beta1*Bt_/Jt;

A = [-Rw/Lw, -Km/Lw, 0, 0;
     Km/Jt, b1, b2, b3;
     0, Kt_*(beta3-1), Kt_*beta2, Kt_*beta1;
     0, beta3, beta2, beta1];

B = transpose([1/Lw, 0, 0, 0]);

C = [Km, 0, 0, 0;
     0, Ns/s, 0, 0;
     0, Ns*Nf/s, 0, 0;
     0, (Bt_+Kt_/s)/Ns/Nf/Nt, 0, (-Kl_-Bl_*s-Bt_*s-Kt_)/(Ns*Nf*Nt)];

D = transpose([0,0,0,0]);

I = eye(size(A));

phi = inv((s*I - A));

M = C*phi*B + D;

testTF = V(1)*M(1)/s;
vt_ = V(2)*M(1);

Q7.Gtau = M(1);

Q7.Gtau = Km*Q4.Ye/(1+(Km)^2*Q4.Ye*Q5.Ym);



Q8.Gq = Km*Q4.Ye*Q5.Ym/(1+Km^2*Q4.Ye*Q5.Ym)/s*180/pi;

%Q8.Gq = testTF*180/pi;




Mtt = 3*Jf*Nf^2/Ns^2+(Mm+Ms)/Ns^2+Jr+Js;

Q9.Gf = Km*Q4.Ye*Q5.Ym/(1+Km^2*Q4.Ye*Q5.Ym)*s*Mtt/Ns;
Q9.Gf = Ns*Q7.Gtau/(3*Nf*Nt);


Lt = Ns^2/(Kt*Nf^2);
Rt = Ns^2/(Bt*Nf^2);
Rl = Ns^2/(Bl*Nf^2*Nt^2);
Ll = Ns^2/(Kl*Nf^2*Nt^2);

zt = (RR(Lt*s,Rt)+RR(Ll*s,Rl));
Q9.Gf = Ns/(Nf*Nt)*Q7.Gtau*Q5.Ym/zt;

Gwr = Q7.Gtau*Q5.Ym;
Gtauf = Ns*Q7.Gtau/(3*Nf);
Gqf = Nf*Gwr/(s*Ns);


Q10.iw = dcgain(Q4.Ye);
Q10.taur = dcgain(Q7.Gtau);
Q10.fs = dcgain(Ns*Q7.Gtau);
Q10.tauf = dcgain(Gtauf);
Q10.ft = dcgain(Ns*Q7.Gtau/(3*Nf*Nt));
Q10.Ktl = dcgain(Gtauf/Gqf);
Q10.qf = dcgain(Gqf)*180/pi;
Q10.ds = dcgain(Gwr/(s*Ns));
Q10.qr = dcgain(Gwr/s)*180/pi;
Q10.vs = Q3.Ks*Q10.qr;


run p1Submit.p


function eq = RR(z1, z2)
    eq = (z1 * z2 ) / (z1 + z2);
end



























%-----------------------Functions 

function [f] = toTransfer(FV, Z, Wn)
    num = FV*Wn^2;
    den = [1 2*Z*Wn Wn^2];
    f = tf(num,den);
end

