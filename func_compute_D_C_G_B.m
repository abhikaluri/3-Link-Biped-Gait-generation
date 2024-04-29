function [D,C,G,B]= func_compute_D_C_G_B(q,dq,param)
%%%%%%  func_compute_D_C_G_B.m
%%%%  04/09/24
%%%%
%%%%
%%%%
%Inputs
q1=q(1);
q2=q(2);
q3=q(3);
%%%%
%%%%
dq1=dq(1);
dq2=dq(2);
dq3=dq(3);
%%%%
%%%%
r=param(1);
m=param(2);
Mh=param(3);
Mt=param(4);
l=param(5);
g=param(6);
%%%%
%%%%
%%%%
%%%%
D=zeros(3,3);
D(1,1) = Mt*l^2 + Mh*r^2 + Mt*r^2 + (5*m*r^2)/4 + (m*r^4)/4 - m*r^2*cos(q2) - 2*Mt*l*r*cos(q3);
D(1,2) = -(m*r^2*(2*cos(q2) - 1))/4;
D(1,3) = Mt*l*(l - r*cos(q3));
D(2,1) = -(m*r^2*(2*cos(q2) - 1))/4;
D(2,2) = (m*r^2)/4;
D(2,3) = 0;
D(3,1) = Mt*l*(l - r*cos(q3));
D(3,2) = 0;
D(3,3) = Mt*l^2;
%%%%
%%%%
C=zeros(3,3);
C(1,1) = (dq2*m*r^2*sin(q2))/2 + Mt*dq3*l*r*sin(q3);
C(1,2) = (m*r^2*sin(q2)*(dq1 + dq2))/2;
C(1,3) = Mt*l*r*sin(q3)*(dq1 + dq3);
C(2,1) = -(dq1*m*r^2*sin(q2))/2;
C(2,2) = 0;
C(2,3) = 0;
C(3,1) = -Mt*dq1*l*r*sin(q3);
C(3,2) = 0;
C(3,3) = 0;
%%%%
%%%%
G=zeros(3,1);
G(1,1) = Mt*g*(l*sin(q1 + q3) - r*sin(q1)) + g*m*((r*sin(q1 + q2))/2 - r*sin(q1)) - Mh*g*r*sin(q1) - (g*m*r^2*sin(q1))/2;
G(2,1) = (g*m*r*sin(q1 + q2))/2;
G(3,1) = Mt*g*l*sin(q1 + q3);
%%%%
%%%%
B=zeros(3,2);
B(1,1) = 0;
B(1,2) = 0;
B(2,1) = 1;
B(2,2) = 0;
B(3,1) = 0;
B(3,2) = 1;
%%%%
%%%%
%%End of code