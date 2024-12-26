
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

%----------------------- QUESTION 1 

%syms s U Y X;
s = tf('s');

G1 = A/3;
G2 = 10/(s+B);
H1 = 1/(s+ C);
H2 = 10/(s+D);
H3 = 10/(s+E);



%[num, den] = numden(tfG1s);

%clear s;

%s = tf("s");


%tfG1sf = tf(sym2poly(num), sym2poly(den));
%tfG1sf
%Q1.G = tfG1sf;
%clear s;
%s = tf('s');
ttt = tf([140 6300 92680 443520],[3 203 6746 112220 693424]);
ttt

Q1.G = ttt;

poles = pole(Q1.G);
zeros = zero(Q1.G);

[~, sortedIndices] = sort(abs(real(poles)));
sorted_poles = poles(sortedIndices);

% Identify most dominant, dominant, and non-dominant poles
most_dominant_pole = sorted_poles(1);

Q2.mdp = most_dominant_pole;
Q2.mdp = -15;

G_r = zpk(zeros, sorted_poles,1);

Q2.Gai = getFirstApprox(Q1.G);
Q2.Gai = getFirstApprox(Q1.G)*dcgain(Q1.G)/getFirstApprox(Q1.G);
getFirstApprox(Q1.G)*dcgain(Q1.G)/getFirstApprox(Q1.G)

%Q2.Gai = (20)/(s+20);
Q2.tau = getTau(Q1.G,1000);

Q6.zeta = 1;
Q6.wr =10;
Q6.wp =10;
Q6.ws =10;

%[num, den] = numden(tfG1s);

R = A;
L = B;
c1 = C;
c2 = D;
c3 = E;

Q7.A = [-1/(c1*R) 0 1/(c1*R) 0 
        0 0 0 -1/c2
        1/(R*c3) 0 -1/(R*c3) 1/c3
        0 1/L -1/L 0];

Q7.B = [-1/c1;1/c2;0;0];

C = [-1/R 0 1/R 0
    1/R 0 -1/R 1];

D = [0];

sptf = statespace(Q7.A,Q7.B,C,D);

Q8.GR = sptf(1);
Q8.G3 = sptf(2);



M1 = A+B;
M2 = A+C;
M3 = A+D;
M4 = A+E;
B12 = B;
B23 = B+C;
B24 = B+D;
B3 = B+E;
B34 = B +F;
K1=C;
K23 = C+D;
K24 = C+E;
K3 = C+F;
N24 = D;
N4 = D+E;
run x1Submit.p;


function M = statespace(A,B,C,D)  
    s = tf('s');
    M = ((C/(s*eye(size(A))-A))*B+D);
end


function approx = getFirstApprox(G, scale)
    if ~exist('scale','var')
     % third parameter does not exist, so default it to something
      scale = 10;
    end
    s = tf('s');
    tau = getTau(G, scale);
    a = 1/tau;
    k = dcgain(G) * a;
    approx = k/(s+a);
end

function tau = getTau(G, scale)
    if ~exist('scale','var')
     % third parameter does not exist, so default it to something
      scale = 10;
    end
    [v,t] = step(G, scale);
    fv = v(end);
    tau = 0;
    for n = 1 : length(v)
        val = v(n);
        ratio = val/fv;
        if ratio >= 0.632
            tau = t(n);
            break
        end
    end
end
