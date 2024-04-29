function [vMh,vMt,vm1,vm2,vcm]= func_compute_vMh_vMt_vm1_vm2_vcm(q,dq,param)
%%%%%%  func_compute_vMh_vMt_vm1_vm2_vcm.m
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
vMh=zeros(2,1);
vMh(1,1) = -dq1*r*cos(q1);
vMh(2,1) = -dq1*r*sin(q1);
%%%%
%%%%
vMt=zeros(2,1);
vMt(1,1) = dq1*(l*cos(q1 + q3) - r*cos(q1)) + dq3*l*cos(q1 + q3);
vMt(2,1) = dq1*(l*sin(q1 + q3) - r*sin(q1)) + dq3*l*sin(q1 + q3);
%%%%
%%%%
vm1=zeros(2,1);
vm1(1,1) = -(dq1*r^2*cos(q1))/2;
vm1(2,1) = -(dq1*r^2*sin(q1))/2;
%%%%
%%%%
vm2=zeros(2,1);
vm2(1,1) = dq1*((r*cos(q1 + q2))/2 - r*cos(q1)) + (dq2*r*cos(q1 + q2))/2;
vm2(2,1) = dq1*((r*sin(q1 + q2))/2 - r*sin(q1)) + (dq2*r*sin(q1 + q2))/2;
%%%%
%%%%
vcm=zeros(2,1);
vcm(1,1) = (dq1*(Mt*(l*cos(q1 + q3) - r*cos(q1)) + m*((r*cos(q1 + q2))/2 - r*cos(q1)) - (m*r^2*cos(q1))/2 - Mh*r*cos(q1)))/(Mh + Mt + 2*m) + (Mt*dq3*l*cos(q1 + q3))/(Mh + Mt + 2*m) + (dq2*m*r*cos(q1 + q2))/(2*(Mh + Mt + 2*m));
vcm(2,1) = (dq1*(Mt*(l*sin(q1 + q3) - r*sin(q1)) + m*((r*sin(q1 + q2))/2 - r*sin(q1)) - Mh*r*sin(q1) - (m*r^2*sin(q1))/2))/(Mh + Mt + 2*m) + (Mt*dq3*l*sin(q1 + q3))/(Mh + Mt + 2*m) + (dq2*m*r*sin(q1 + q2))/(2*(Mh + Mt + 2*m));
%%%%
%%%%
%%End of code