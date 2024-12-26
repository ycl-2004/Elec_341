%%%%%%%%%%%%%%%%%%%%
% Assignment 2 ELEC 341
%
% Author: YI-CHEN LIN

% FV -> Final Value
% OS -> Overshoot (Peak-FV)/FV *100%
% Tr -> Rise_Time, time it takes from 10% to 90% of the final value
% Tp -> Peak_Time, time it takes to reach the first maximum overshoot
% Ts -> Settle_Time, time it takes for system to stable

%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

SN = 48861785;
a2DSPlot(48861785)

%--------------------------------------------------------

Q1.FV =51.5; %49 50 1.5/3    51 52 -> 3/3 
z1f = Q1.FV;

overshoot1 = (80.2-Q1.FV)/Q1.FV; %78.8 -> 1.5/3, 80~80.2 ->3,
Q1.OSy =overshoot1*100;

Q1.zeta =sqrt(log(overshoot1)^2/(pi^2+log(overshoot1))); %Formula Sheet

Q1.Tr =0.027; % 0.022 to 0.032 <20% G,
Q1.Tp =0.047; % 0.042 to 0.052 <20% G

beta1 = sqrt(1-Q1.zeta^2);

Q1.wn =pi/(Q1.Tp*beta1);

Q1.Ts =4/(Q1.zeta*Q1.wn);

%----------------------------- Question 1 End 
FV2 = 33.55; % 33.3~.8 error <3%
overshoot2 = (49.6264-FV2)/FV2;

Q2.OSy = overshoot2*100;

Q2.Tr = 0.144; %  0.13<to<0.158 G


Q2.zeta = -log(overshoot2) / sqrt(pi^2 + (log(overshoot2))^2); %Similar to previous one just w negative sign

beta2 = sqrt(1-Q2.zeta^2);

Q2.wnr = (pi-atan(beta2/Q2.zeta))/(Q2.Tr*beta2);


s = tf('s');
xfer2 = toTransfer(FV2,Q2.zeta,Q2.wnr);

Q2.G = xfer2;


%----------------------------- Question 2 End

Q3.Tp = 0.373; % 0.355<to<0.391 G

Q3.wnp = pi/(Q3.Tp*beta2);

xfer3 = toTransfer(FV2,Q2.zeta,Q3.wnp);

Q3.G = xfer3;

%----------------------------- Question 3 End

Q4.Ts =0.748; % 0.711 0.785 G
Q4.wns = (3.9+log(1/beta2))/(Q4.Ts*Q2.zeta);

xfer4 = toTransfer(FV2,Q2.zeta,Q4.wns);
Q4.G =xfer4 ;

%----------------------------- Question 4 End

settle_time_5 = 0.705; % 0.64<to<0.77 G
Q5.tau =settle_time_5/4; 

Q5.Tr1 = 0.28;

wn5 = pi/0.63; % 0.57<to<0.69 G where 0.6 -> Tp

Q5.zeta = 1/(wn5*Q5.tau);


%----------------------------- Question 5 End

Q6.Tr1 = 0.31; % 0.26<to<0.31 G
Q6.wn1 = (4.44*Q5.zeta-1.15)/Q6.Tr1;

xfer6 = toTransfer(540,Q5.zeta,Q6.wn1);
Q6.G = xfer6;

%----------------------------- Question 6 End

%---------------------Function Used

function [f] = toTransfer(FV, ZETA, Wn)
    num = FV*Wn^2;
    den = [1 2*ZETA*Wn Wn^2];
    f = tf(num,den);
end