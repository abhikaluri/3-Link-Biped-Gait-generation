function [beta2]= func_compute_beta2(q,dq,param)
%%%%%%  func_compute_beta2.m
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
beta2=zeros(1,1);
beta2(1,1) = (m*(2*(dq1*((r*cos(q1 + q2))/2 - r*cos(q1)) + (dq2*r*cos(q1 + q2))/2)*((r*cos(q1 + q2))/2 - r*cos(q1)) + 2*((r*sin(q1 + q2))/2 - r*sin(q1))*(dq1*((r*sin(q1 + q2))/2 - r*sin(q1)) + (dq2*r*sin(q1 + q2))/2)))/2 + (Mh*(2*dq1*r^2*cos(q1)^2 + 2*dq1*r^2*sin(q1)^2))/2 + (m*((dq1*r^4*cos(q1)^2)/2 + (dq1*r^4*sin(q1)^2)/2))/2 + (Mt*(2*(dq1*(l*cos(q1 + q3) - r*cos(q1)) + dq3*l*cos(q1 + q3))*(l*cos(q1 + q3) - r*cos(q1)) + 2*(l*sin(q1 + q3) - r*sin(q1))*(dq1*(l*sin(q1 + q3) - r*sin(q1)) + dq3*l*sin(q1 + q3))))/2;
%%%%
%%%%
%%End of code