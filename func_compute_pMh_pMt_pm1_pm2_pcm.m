function [pMh,pMt,pm1,pm2,pcm]= func_compute_pMh_pMt_pm1_pm2_pcm(q,dq,param)
%%%%%%  func_compute_pMh_pMt_pm1_pm2_pcm.m
%%%%  02/25/24
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
pMh=zeros(2,1);
pMh(1,1) = r*sin(q1);
pMh(2,1) = r*cos(q1);
%%%%
%%%%
pMt=zeros(2,1);
pMt(1,1) = r*sin(q1) - l*sin(q1 - q3);
pMt(2,1) = r*cos(q1) - l*cos(q1 - q3);
%%%%
%%%%
pm1=zeros(2,1);
pm1(1,1) = (r*sin(q1))/2;
pm1(2,1) = (r*cos(q1))/2;
%%%%
%%%%
pm2=zeros(2,1);
pm2(1,1) = r*sin(q1) - (r*sin(q1 - q2))/2;
pm2(2,1) = r*cos(q1) - (r*cos(q1 - q2))/2;
%%%%
%%%%
pcm=zeros(2,1);
pcm(1,1) = (mt*(r*sin(q1) - l*sin(q1 - q3)) + m*(r*sin(q1) - (r*sin(q1 - q2))/2) + (m*r*sin(q1))/2 + mh*r*sin(q1))/(2*m + mh + mt);
pcm(2,1) = (mt*(r*cos(q1) - l*cos(q1 - q3)) + m*(r*cos(q1) - (r*cos(q1 - q2))/2) + (m*r*cos(q1))/2 + mh*r*cos(q1))/(2*m + mh + mt);
%%%%
%%%%
%%End of code