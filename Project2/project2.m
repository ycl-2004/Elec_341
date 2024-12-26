%%%%%%%%%%%%%%%%%%%%
% Project_1 ELEC 341
%
% Author: YI-CHEN LIN

%%%%%%%%%%%%%%%%%%%%
clear all;
clc

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

SN = 48861785;

SN = 48861785;


SN = 48861785;
%p1DSPlot(SN);
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

tp = 0.001705;
fv = 18;
peak = 22.5;

%Q1.Ga = toTransfer(fv,zeta,wnp);

%-----------------------Question 2 Underdamp
zeta2 = -log(overshoot1) / sqrt(pi^2 + (log(overshoot1))^2); 
beta2 = sqrt(1-zeta2^2);

wn2 = (pi-atan(beta2/zeta2))/(Q1.Tr*beta2);

os = (peak-fv)/fv;
zeta = sqrt(log(os)^2/(pi^2+log(os)^2));
beta = sqrt(1-zeta^2);
wnp = pi/(tp*beta);

Q2.Ga = toTransfer(fv,zeta,wnp);

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



%------------------------Question 4
s = tf("s");

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

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

%------------------------ PROJECT 2 START--------------------------

%Begin with 0 order controller, Gc = 1;

A = 4 + 10;
B = 8 + 10;
C = 8 + 10;
D = 6 + 10;
E = 1 + 10;
F = 7 + 10;
G = 8 + 10;
H = 5 + 10;

nn = Ns/(Nt*Nf);
Q11.G = Q2.Ga*Q8.Gq*pi/180;
Q11.Ktj = nn;  %feedback(Q11.G,Hs Hc);   %Ktj = x1/yd
Q11.Kjt = 1/Q11.Ktj;   %Kjt = ya/x4 


CF = A*B*5;
DC = (C + D + E + F)/100;
deltat = 1/CF;

Q12.Noh = DC + 0.5;


poles = pole(Q11.G*Q3.Ds);
format long;
%disp(poles);

Q13.wd = 0.95; %1 2<4% 
Q13.Nf = CF/(Q13.wd*10+DC);
Q13.tau = (Q13.Nf+0.5)/CF;
Q13.beta = exp(-deltat/Q13.tau);
Q13.N = Q13.Nf + 0.5;

Nterm = 4*Q13.tau/deltat+1;
Q14.NC = 501;

t_values =1:1:Q14.NC;
Craw_values = zeros(size(t_values));
total = 0;

for i = 1:length(t_values)
     %Craw_values(i) = exp(-t_values(i)/tauf);
     Craw_values(i) = exp(-4*t_values(i)/(Q14.NC-1.7));
     total = total+Craw_values(i);
end

Q14.W = Craw_values/total;

Hs = Q3.Ds*Q3.Ks;
Q15.Hc =1/(Q13.tau*s+1)/Q3.Ks*pi/180; 
Q15.H = Q15.Hc*Hs*180/pi;

Dnew = CF/((Q13.N)*s/2+CF)/s;
Q16.Dp = Dnew;
[Q16.K0,pm,wx,wcp] = margin(Q16.Dp*Q11.G*Q15.H);
[K0,pm,Q16.wxo,wcp] = margin(Q16.K0*Q16.Dp*Q11.G*Q15.H);

[zz1,zz2,pm] = optimizeDoubleZero(-Q16.wxo,Q16.K0*Q16.Dp*Q11.G*Q15.H,0.1, deg2rad(1));

Q17.Z = [zz1,zz2];
Q17.PM = pm;

Dz = (s-zz1)*(s-zz2)/(zz1*zz2);
Q17.D = Dz*Q16.Dp;

K_val = 0.284
[K0,pm,wxo,wcp] = margin(K_val*Q17.D*Q11.G*Q15.H);
Q18.K = K_val;

p = -1*CF/(Q13.N/2);

Q19.Kp = K_val*(1/p-(zz1+zz2)/(zz1*zz2));
Q19.Ki = K_val;
p = -1*CF*Q13.N;
Q19.Kd = K_val*(1/p^2-(zz1+zz2-p)/(p*zz1*zz2));
Q19.Kd = 0.035;

K = 1;

PID = K*(Q19.Kp + Q19.Ki*1/s + 0.5*Q19.Kd*(-p*s)/(s-p));

% Kp + , Ki + ,Kd -
% ki -> steady large kp->t=0

%a = 1;
%b = 1;
%c = 0.5;
%PID = K*(a*Q19.Kp + b*Q19.Ki*1/s + c*Q19.Kd*(-p*s)/(s-p));

%PID = Q18.K * (Q19.Kp + Q19.Ki*1/s + Q19.Kd*(2*CF*s)/(Q13.N*s+2*CF))
Gcl = feedback(Q18.K*Q17.D*Q11.G, Q15.H);

%step(Gcl);

Q20.Tr = 0.17;
Q20.Tp = 0.27;
Q20.Ts = 3.23; %2.4225
Q20.OSu = 25;  %15
Q20.OSy = 25;
Q20.Ess = 0;

K_tune = 0.99; %0.284 NEXT TIME TRY 0.6 up by 0.01
Gcl = feedback(K_tune*PID*Q11.G, Q15.H);
step(Gcl)
hold on

yline(1.02)
yline(0.98)

Q21.K = K_tune;
Q21.Kp = Q19.Kp;
Q21.Ki = Q19.Ki;
Q21.Kd = Q19.Kd*0.5;

Q22.Tr = 0.35; % 75~150 0.6
Q22.Tp = 0.41; % 0.38 ~ 10~20E
Q22.Ts = 1.7; % 10~25 1.9
Q22.OSu = 12;
Q22.OSy = 12;
Q22.Ess = 0;

run p2Submit.p



function eq = RR(z1, z2)
    eq = (z1 * z2 ) / (z1 + z2);
end



function [zero1,zero2,PM] = optimizeDoubleZero(start,pdy,resR, resA)
% optimizeDoubleZero Optimizes two zeros for given partial dynamics system.
% Author: Henrique Saito
% 
%   PARAMS:
%
%       start: Real number representing initial point of search. Should be
%       negative of crossover frequency.
%       
%       pdy: Partial dynamics LTI system. pdy = K0*Dp*Gp*H
%       
%       resR (def: 0.1): Search resolution for real values and polar
%       radius.
%
%       resA (def: deg2rad(1)): Search resolution for polar angle. In
%       radians.

    if ~exist('resR','var')
      resR = 0.1; % default resolution
    end
    if ~exist('resA', 'var')
        resA = deg2Rad(1);
    end

    [dzero, ~] = optimizeRealDoubleZero(start,pdy,resR);
    [sep, ~, ~] = optimizeSeparation(dzero, 0, pdy, resR);
    if sep > 0
        [zero1,zero2,PM] = optimizeRealOnly(dzero, sep, pdy, resR);
    else
        [zero1,zero2,PM] = optimizeComplex(-1*dzero,0,pdy,resR, resA);
    end % if
end


function [z,pm] = optimizeRealDoubleZero(start,pdy,res)
    currPM = getPM(start,start,pdy);
    zero = start;
    done = false;
    haveChecked = [zero];
    while ~done
        inc = zero + res;
        if ~isin(haveChecked, inc)
            haveChecked(end+1) = inc;
            incPM = getPM(inc,inc,pdy);
            if (incPM > currPM)
                zero = inc;
                currPM = incPM;
                continue
            end % if
        end % if
        dec = zero - res;
        if ~isin(haveChecked, dec)
            haveChecked(end+1) = dec;
            decPM = getPM(dec,dec,pdy);
            if (decPM > currPM)
                zero = dec;
                currPM = decPM;
                continue
            end % if
        end % if
        done = true;
    end % while
    z = zero;
    pm = currPM;
    disp(['Double real zero: ', num2str(zero), ' | PM: ', num2str(pm)])
end

function [separation,pm,changed] = optimizeSeparation(center, initSep, pdy, res)
    z1 = center + initSep;
    z2 = center - initSep;
    sep = initSep;
    currPM = getPM(z1,z2,pdy);
    done = false;
    changed = false;
    haveChecked = [initSep];
    while ~done
        incSep = sep + res;
        if ~isin(haveChecked, incSep)
            haveChecked(end+1) = incSep;
            incZ1 = center + incSep;
            incZ2 = center - incSep;
            if (incZ1 <= 0)
                incPM = getPM(incZ1,incZ2,pdy);
                if (incPM > currPM)
                    changed = true;
                    sep = incSep;
                    currPM = incPM;
                    continue
                end % if
            end % if
        end % if
        decSep = sep - res;
        if ~isin(haveChecked, decSep)
            if (decSep < 0)
                decSep = 0;
            end % if
            haveChecked(end+1) = decSep;
            decZ1 = center + decSep;
            decZ2 = center - decSep;
            decPM = getPM(decZ1, decZ2, pdy);
            if (decPM > currPM)
                changed = true;
                sep = decSep;
                currPM = decPM;
                if (sep == 0)
                    done = true;
                end % if
                continue
            end % if
        end % if
        done = true;
    end % while
    separation = sep;
    pm = currPM;
    disp(['Separation: ',num2str(sep), ' | PM: ', num2str(pm)]);
end
    
function [zero1,zero2,PM] = optimizeRealOnly(initCenter, initSep, pdy, res)
    center = initCenter;
    sep = initSep;
    changed = true;
    tuningCenter = true;
    while changed
        if tuningCenter
            [center, pm, changed] = optimizeRealCenter(center, sep, pdy, res);
            tuningCenter = false;
        else
            [sep, pm, changed] = optimizeSeparation(center, sep, pdy, res);
            tuningCenter = true;
        end % if
    end % while
    zero1 = center + sep;
    zero2 = center - sep;
    PM = pm;
    disp(['Optimized Real Zeros: ', num2str(zero1), ' ', num2str(zero2), ' | PM: ', num2str(pm)]);
end

function [centerOut, pm, changed] = optimizeRealCenter(initCenter, sep, pdy, res)
    z1 = initCenter + sep;
    z2 = initCenter - sep;
    center = initCenter;
    currPM = getPM(z1,z2,pdy);
    done = false;
    changed = false;
    haveChecked = [initCenter];
    while ~done
        incCenter = center + res;
        if ~isin(haveChecked, incCenter)
            haveChecked(end+1) = incCenter;
            incZ1 = incCenter + sep;
            incZ2 = incCenter - sep;
            if (incZ1 <= 0)
                incPM = getPM(incZ1,incZ2,pdy);
                if (incPM > currPM)
                    changed = true;
                    center = incCenter;
                    currPM = incPM;
                    continue
                end % if
            end
        end % if
        decCenter = center - res;
        if ~isin(haveChecked, decCenter)
            haveChecked(end+1) = decCenter;
            decZ1 = sep + decCenter;
            decZ2 = decCenter - sep;
            decPM = getPM(decZ1, decZ2, pdy);
            if (decPM > currPM)
                changed = true;
                center = decCenter;
                currPM = decPM;
                continue
            end % if
        end % if
        done = true;
    end % while
    centerOut = center;
    pm = currPM;
    disp(['Center: ',num2str(center), ' | PM: ', num2str(pm)]);  
end

function [angle, pm, changed] = optimizeComplexAngle(radius, initAngle, pdy, res)
    ang = initAngle;
    ang1 = pi - ang;
    z1 = radius.*exp(ang1*1i);
    z2 = conj(z1);
    currPM = getPM(z1,z2,pdy);
    done = false;
    changed = false;
    haveChecked = [initAngle];
    while ~done
        incAng = ang + res;
        if ~isin(haveChecked, incAng)
            haveChecked(end+1) = incAng;
            incAng1 = pi - incAng;
            incZ1 = radius.*exp(incAng1*1i);
            incZ2 = conj(incZ1);
            incPM = getPM(incZ1,incZ2,pdy);
            if (incPM > currPM)
                changed = true;
                ang = incAng;
                currPM = incPM;
                continue
            end % if
        end % if
        decAng = ang - res;
        if ~isin(haveChecked, decAng)
            if (decAng < 0)
                decAng = 0;
            end % if
            haveChecked(end+1) = decAng;
            decAng1 = pi - decAng;
            decZ1 = radius.*exp(decAng1*1i);
            decZ2 = conj(decZ1);
            decPM = getPM(decZ1, decZ2, pdy);
            if (decPM > currPM)
                changed = true;
                ang = decAng;
                currPM = decPM;
                if (ang == 0)
                    done = true;
                end % if
                continue
            end % if
        end % if
        done = true;
    end % while
    angle = ang;
    pm = currPM;
    disp(['Angle: ', num2str(ang), ' | PM: ',  num2str(pm)]);
end

function [radius, pm, changed] = optimizeComplexRadius(initRadius, angle, pdy, res)
    rad = initRadius;
    ang1 = pi - angle;
    z1 = rad.*exp(ang1*1i);
    z2 = conj(z1);
    currPM = getPM(z1,z2,pdy);
    done = false;
    changed = false;
    haveChecked = [initRadius];
    while ~done
        incRad = rad + res;
        if ~isin(haveChecked, incRad)
            haveChecked(end+1) = incRad;
            incZ1 = incRad.*exp(ang1*1i);
            incZ2 = conj(incZ1);
            incPM = getPM(incZ1,incZ2,pdy);
            if (incPM > currPM)
                changed = true;
                rad = incRad;
                currPM = incPM;
                continue
            end % if
        end % if
        decRad = rad - res;
        if ~isin(haveChecked, decRad)
            if (decRad < 0)
                decRad = 0;
            end % if
            haveChecked(end+1) = decRad;
            decZ1 = decRad.*exp(ang1*1i);
            decZ2 = conj(decZ1);
            decPM = getPM(decZ1, decZ2, pdy);
            if (decPM > currPM)
                changed = true;
                rad = decRad;
                currPM = decPM;
                if (rad == 0)
                    done = true;
                end % if
                continue
            end % if
        end % if
        done = true;
    end % while
    radius = rad;
    pm = currPM;
    disp(['Radius: ',num2str(radius), ' | PM: ', num2str(pm)]);  
end

function [zero1, zero2, PM] = optimizeComplex(initRad, initAngle, pdy, resR, resA)
    rad = initRad;
    ang = initAngle;
    changed = true;
    tuningAngle = true;
    first = true;
    while changed | first
        if tuningAngle
            [ang, pm, changed] = optimizeComplexAngle(rad, ang, pdy, resA);
            tuningAngle = false;
        else
            [rad, pm, changed] = optimizeComplexRadius(rad, ang, pdy, resR);
            tuningAngle = true;
            first = false;
        end % if
    end % while
    ang1 = pi - ang;
    zero1 = rad.*exp(ang1*1i);
    [zero1, changed] = optimizeCartesian(zero1, pdy, resR);
    if changed
        [zero1,zero2,PM] = optimizeComplex(abs(zero1), pi-angle(zero1), pdy, resR, resA);
    else
        zero2 = conj(zero1);
        PM = pm;
        disp(['Optimized Complex Zeros: ', num2str(zero1), ' ', num2str(zero2), ' | PM: ', num2str(pm)]);
    end
end

function [zero1,changed] = optimizeCartesian(init,pdy,res)
    disp(['Beginning Cartesian Tune: ', num2str(init)])
    zero = init;
    tuningReal = true;
    changed = true;
    changed2 = false;
    first = true;
    while changed | first
        if tuningReal
            [zero,changed] = optimizeCartesianReal(zero, pdy, res);
            tuningReal = false;
            if changed
                changed2 = true;
            end
        else
            [zero,changed] = optimizeCartesianImaginary(zero, pdy, res);
            tuningReal = true;
            first = false;
        end
    end
    zero1 = zero;
    zero2 = conj(zero);
    changed = changed2;
    disp(['Optimized Cartesian Complex Zeros: ', num2str(zero1), ' ', num2str(zero2)]);
end

function [out_zero, changed] = optimizeCartesianReal(zero_init, pdy, res)
    currPM = getPM(zero_init,conj(zero_init),pdy);
    zero = zero_init;
    done = false;
    haveChecked = [real(zero)];
    while ~done
        inc = zero + res;
        if ~isin(haveChecked, real(inc))
            haveChecked(end+1) = real(inc);
            incPM = getPM(inc,conj(inc),pdy);
            if (incPM > currPM)
                zero = inc;
                currPM = incPM;
                continue
            end % if
        end % if
        dec = zero - res;
        if ~isin(haveChecked, real(dec))
            haveChecked(end+1) = real(dec);
            decPM = getPM(dec,conj(dec),pdy);
            if (decPM > currPM)
                zero = dec;
                currPM = decPM;
                continue
            end % if
        end % if
        done = true;
    end % while
    out_zero = zero;
    pm = currPM;
    changed = out_zero ~= zero_init;
    disp(['Cartesian real zero: ', num2str(zero), ' | PM: ', num2str(pm)])
end

function [out_zero, changed] = optimizeCartesianImaginary(zero_init, pdy, res)
    currPM = getPM(zero_init,conj(zero_init),pdy);
    zero = zero_init;
    done = false;
    haveChecked = [imag(zero)];
    while ~done
        inc = zero + i*res;
        if ~isin(haveChecked, imag(inc))
            haveChecked(end+1) = imag(inc);
            incPM = getPM(inc,conj(inc),pdy);
            if (incPM > currPM)
                zero = inc;
                currPM = incPM;
                continue
            end % if
        end % if
        dec = zero - i*res;
        if ~isin(haveChecked, imag(dec))
            haveChecked(end+1) = imag(dec);
            decPM = getPM(dec,conj(dec),pdy);
            if (decPM > currPM)
                zero = dec;
                currPM = decPM;
                continue
            end % if
        end % if
        done = true;
    end % while
    out_zero = zero;
    pm = currPM;
    changed = out_zero ~= zero_init;
    disp(['Cartesian imaginary zero: ', num2str(zero), ' | PM: ', num2str(pm)])
end

function pm = getPM(z1,z2,pdy)
    s = tf('s');
    Dz = (s-z1)*(s-z2)/(z1*z2);
    openloop = Dz * pdy;
    [bbb,~] = phsMargin(openloop);
    pm = abs(bbb);
end

function [pm wcp] = phsMargin(KDGH)

  % Find pm x-over freq
  [x pm x wcp] = margin(KDGH);            % x is a dummy variable

  % Use bode to re-calculate pm WRT -180 (not -540)
  if ~isnan(wcp)
    [x phase] = bode(KDGH, wcp);          % x is a dummy variable
    pm = phase + 180;
  end

end % function


function out = isin(array, element)
    arr_ismember = ismembertol(array,element);
    out = sum(arr_ismember) > 0;
end
















function Gp = getPlantAdmittances(Ym,Ye,Km, Ga)
    if ~exist('Ga','var')
        Ga = 1;
    end
    Gp = Ga*feedback(Ye*Km*Ym, Km);
end



%-----------------------Functions 

function [f] = toTransfer(FV, Z, Wn)
    num = FV*Wn^2;
    den = [1 2*Z*Wn Wn^2];
    f = tf(num,den);
end

