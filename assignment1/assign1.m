%%%%%%%%%%%%%%%%%%%%
% Intro
%
% Author: YICHEN LIN
%%%%%%%%%%%%%%%%%%%%
clear all; 
clc;

SN = 48861785;

% PaR1meters
R1 = 140; % (Ohm)
R2 = 160; % (Ohm)
L1 = 18e-3; % (H)
L2 = 11e-3; % (H)
C1 = 18e-6; % (F)
C2 = 17e-6; % (F)

s = tf('s'); 

Q1.G = useTF(R1, R2, L1,L2,C1, C2); 
Q1.G
%------------------------------For Question 2

s = tf('s'); 

zr1 = R1;
zr2 = R2;

zl1= (L1*s);
zl2 =(L2*s);

zc1 = 1/(s*C1);
zc2 = 1/(s*C2);

zt1 = zr2 +zl2+zc2;
zt4 = RR(zr1+zc1,zl1);

yt1 = 1/zt1;
yt2 = 1/zr2;
yt3 = 1/zl2;
yt4 = (s^2*L1*C1+s*C1*R1+1)/(s*L1*(1+s*R1*C1));
yt5 = 1/zc2;

Y_matrix = [yt1+yt2 -yt2 -yt1
            -yt2 yt2+yt3+yt4 -yt3
            -yt1 -yt3 yt1+yt3+yt5];

I_vector = [1;0;0];
V_vector = inv(Y_matrix)*I_vector;



Q2.G = V_vector(2)*yt4;

%%useTF2(R1, R2, L1,L2,C1, C2); 
%%Q2.G

%------------------------------For Question 3

Q3.Kdc = 61;
Q3.sigma = -22;

%------------------------------Q3, Envolope Calculation

d_factor = -22;
k = Q3.Kdc;
t = linspace(0, 0.3,100);

numer = 150^2 * k;
denom = [1 -2*d_factor 150^2];
G = tf(numer, denom);

pos_env =  k*(1+exp(d_factor .* t));
neg_env =  k*(1-exp(d_factor .* t));

t_step = linspace(0, 0.3, 1000); % Time vector for step response

% Get step response
[response, t_step] = step(G, t_step);

% Convert step response time to milliseconds
t_step = t_step * 1000;

%------------------------ For Step Response (Question 4)

%------------------------ Combined Plot

%figure;
%hold on;
%plot(t * 1000, pos_env, 'r', 'LineWidth', 1.5); 
%plot(t * 1000, neg_env, 'r', 'LineWidth', 1.5);
%plot(t_step, response, 'k', 'LineWidth', 3);
%hold off;

w_n = 150;           
d_factor = -22;    
k = Q3.Kdc;

G = (k * w_n^2) / (s^2 - 2*(d_factor)*s + w_n^2);
%figure;
%step(G);
%grid on;
%hold off;

poles_G = pole(G);

pole_min = real(poles_G(1));
imp_pole = 10*pole_min;

Q4.p = imp_pole;

t = 0:1e-3:0.1;

G_impulse = G*(-1*imp_pole / (s - imp_pole));
[v_a, t_1] = impulse(G_impulse, t);

Q4.y = v_a;

%------------------------------

Q5.tau = 0.005;

%u = (t <= Q5.tau) * (1/Q5.tau);
%u
%U = tf([1 -exp(-Q5.tau)], [Q5.tau 0]);
%U

tf5 = 200*(1/s-exp(-1*s/200)/s);

GG = tf5*G;


% Store the result in Q5.y
Q5.y = impulse(GG,t);

%--------------------------
%-----------------------------------------------------------------

function [num, den] = numDen(R1, R2, L1,L2,C1, C2)
num = [C1*L1*L2*R1 C1*L1*R1*R2 L1*R2];
den = [C1*C2*L1*L2*R1*R2 C1*L1*L2*R1+C2*L1*L2*R2 C2*L2*R1*R2 L2*R1];
end

function xfer = useTF(R1, R2, L1,L2,C1, C2)

[num1, den1] = numDen(R1, R2, L1,L2,C1, C2);
xfer = tf(num1, den1);
end % function

%-------------------------------End Question 1


%-------------------------------Question 2 Function
function xfer = useTF2(R1, R2, L1,L2,C1, C2)

[num2, den2] = numDen2(R1, R2, L1,L2,C1, C2);
xfer = tf(num2, den2);
end % function


function [num2, den2] = numDen2(R1, R2, L1,L2,C1, C2)
num2 = [C1*C2^2*L1*L2^2, ...
    C1*C2^2*L1*L2*R2 + C1*C2^2*L2^2*R1, ...
    C1*C2^2*L2*R1*R2 + 3*C1*C2*L1*L2, ...
    2*C1*C2*L1*R2 + 3*C1*C2*L2*R1 + C2^2*L2*R2, ...
    2*C1*C2*R1*R2 + C1*L1 + C2^2*L2^2 + 3*C2*L2, ...
    C1*R1 + 2*C2*R2, ...
    1];

den2 =  [C1*C2^2*L1*L2^2, ...
    C1*C2^2*L2^2*R1 + 2*C1*C2^2*L1*L2*R1 + 2*C1*C2^2*L1*L2*R2, ...
    2*C1*C2^2*L1*R1*R2 + 2*C1*C2^2*L2*R1*R2 + C2^2*L2^2 + 2*C2^2*L1*L2 + 3*C1*C2*L1*L2, ...
    2*C2^2*L1*R2 + 2*C2^2*L2*R2 + C1*C2*L1*R1 + 2*C1*C2*L1*R2 + 3*C1*C2*L2*R1, ...
    2*C1*C2*R1*R2 + C1*L1 + C2*L1 + 3*C2*L2, ...
    C1*R1 + 2*C2*R2, ...
    1];

end

%-----------------------------Question 2 End

function resul = RR(a,b)
    resul = (a+b)/(a*b);
end

run a1Submit.p;