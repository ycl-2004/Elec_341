M0 = (4+10)/5;
M1 = (8+10)/10;
M2 = (8+10)/10;
M3 = (6+10)/5;

B20 = (1+10)/2;
B21 = (7+10)/3;
B31 = (8+10)/4;

K0 = (4+10);
K1 = (8+10);
K20 = (8+10);
K32 = (6+10)/3; 



Q1.A = [-B20/M0, 0, B20/M0, 0, -1/M0, 0, 1/M0, 0 
     0, (-B21-B31)/M1, B21/M1, B31/M1, 0, -1/M1, 0, 0
     B20/M2, B21/M2, (-B20-B21)/M2, 0,0,0,-1/M2, 1/M2
     0, B31/M3, 0, -B31/M3, 0, 0, 0, -1/M3
     K0, 0, 0, 0, 0, 0, 0, 0
     0, K1, 0, 0, 0, 0, 0, 0
     -K20, 0, K20, 0, 0, 0, 0, 0
     0, 0, -K32, K32, 0, 0, 0, 0];

Q1.B = [1/M0;0;0;0;0;0;0;0];

C = [0 B31/M3 0 -B31/M3 0 0 0 -1/M3
    0 -B21 B21 0 0 0 0 0];

s = tf('s');

C = [0 0 0 1/s 0 0 0 0
    0 B21 -B21 0 0 0 0 0];
D = [0
     0];


Gd3 = ((C/(s*eye(size(Q1.A))-Q1.A))*Q1.B+D);
%Gd3 = tf(ss(Q1.A,Q1.B,C,D));

Q2.Gd3 =Gd3(1);
Q2.Gf21 =-Gd3(2);