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