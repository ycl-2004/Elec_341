% Params
clear all;
clc;

A_sn = 1 + 10;
B_sn = 9 + 10;
C_sn = 1 + 10;
D_sn = 7 + 10;
E_sn = 7 + 10;
F_sn = 0 + 10;
G_sn = 6 + 10;
H_sn = 2 + 10;

Hs_gain = 5/F_sn;
Gp_gain = 3/G_sn;
Hs_pole = -15*E_sn;
Gp_p1 = -3*B_sn;
Gp_p2 = -A_sn;
Gp_p3 = (-1+1i)*2*D_sn;
Gp_p4 = (-1-1i)*2*D_sn;
Gp_zero = -5*C_sn;

Kh = 1/Hs_gain;

Gp_mod = zpk(Gp_zero, [Gp_p1, Gp_p2, Gp_p3, Gp_p4], 1);
Gp_mod_k = dcgain(Gp_mod);
Gp_k_e = Gp_gain/Gp_mod_k;
Gp = Gp_mod * Gp_k_e;

Hs_mod = zpk([], Hs_pole, 1);
Hs_mod_k = dcgain(Hs_mod);
Hs_k_e = Hs_gain/Hs_mod_k;
Hs = Hs_mod * Hs_k_e;

GH = Gp*Hs*Kh

Ku = 36.1;

K = Ku/2;

Gcl = feedback(K*Gp, Hs*Kh);

Q1.Kh = Kh;
Q1.GH = GH;
Q1.Ku = Ku;


Q2.Gcl = Gcl;

step(Q2.Gcl)
Q2.Ts = 0.4840;
Q2.Ess = 100-77.5;
Q2.Pos = (1.1174-0.775)/(0.775)*100;
Q2.Gos = (1.1174-1)*100;

Klnt = Kl*Nt^2;
Blnt = Bl*Nt^2;

Klnf = Klnt*Nf^2;
Blnf = Blnt*Nf^2;
Ktnf = Kt*Nf^2;
Btnf = Bt*Nf^2;
Bfnf = Bf*Nf^2;
Jrnf = Jr*Nf^2;

Bsns = Bs/Ns^2;
Mmtns = (Mm+Ms)/Ns^2;
Jrns = Jrnf/Ns^2;
Bfns = Bfnf/Ns^2;
Ktns = Ktnf/Ns^2;
Btns = Btnf/Ns^2;
Klns = Klnf/Ns^2;
Blns = Blnf/Ns^2;

Jt = Jr+Js;
Br = Br;

Rbsns = 1/Bsns;
Cmmtns = Mmtns;
Cjrns = Jrns;
Rbfns = 1/Bfns;
Lktns = 1/Ktns;
Rbtns = 1/Btns;
Lklns = 1/Klns;
Rblns = 1/Blns;
Cjt = 1/Jt;
Rbr = 1/Br;

zRbsns = Rbsns;
zCmmtns = 1/(Cmmtns*s);
zCjrns = 1/(Cjrns*s);
zRbfns = Rbfns;
zLktns = Lktns*s;
zRbtns = Rbtns;
zLklns = Lklns*s;
zRblns = Rblns;
zCjt = 1/(Cjt*s);
zRbr = Rbr;

zt1n = zRblns*zLklns/(zLklns+zRblns) + zLktns*zRbtns/(zRbtns+zLktns);

zt2n = 1/zt1n + 1/zRbfns + 1/zCjrns + 1/zRbsns + 1/zCmmtns + 1/zRbr;
zt2n = 1/zt2n;

zt3n = 1/(1/zCjt);

izt3n = zt2n/(zt2n+zt3n);
vzt3n = izt3n*zt3n;

zRl = 1/Blns;
zLl = 1/Klns*s;
zRt = 1/Btns;
zLt = 1/Ktns*s;
zRf = 1/Bfns;
zJr = 1/(s*Jrns);
zRs = 1/Bsns;
zMt = 1/(s*Mmtns);
zJt = 1/(s*Jt);
zRr = 1/Br;


zt1 = zRl*zLl/(zRl+zLl) + zLt*zRt/(zLt+zRt);

zt2 = 1/zt1 + 1/zRf + 1/zJr + 1/zRs + 1/zMt + 1/zJt;
zt2 = 1/zt2;

ir2 = 1*zt2/(zt2+zRr);
ww = ir2*zRr;

Q5.Ym = vzt3n;
Q5.Ym = minreal(1/(s*Jtotal + Btotal + Ktotal/s))*(Ns*Nf*Nt);
Q5.Ym