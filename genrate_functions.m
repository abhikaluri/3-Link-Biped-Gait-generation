
% Define symbolic variables
syms q1 q2 q3 p_h p_v dp_h dp_v dq1 dq2 dq3 real
syms r m Mh Mt l g real
% Define parameters in a vector
params = [r,m,Mh,Mt,l,g];
set_path
% Define generalized coordinates
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
% Position vectors
pMh = [r * sin(-q1); r * cos(-q1)];
pMt = pMh + [l * sin(-q1+pi-q3);  l * cos(-q1+pi-q3)];
pm1 = r*0.5 * pMh;
pm2 = pMh + r/2 * [sin(q2 + q1); -cos(q2 + q1)];
pcm = (m*pm1 +m*pm2+Mh*pMh+Mt*pMt)/(m+m+Mh+Mt);
P2 = pMh + [r * sin(q2 + q1); -r * cos(q2 + q1)];

write_symbolic_term_to_mfile(q,dq,params,pMh,pMt,pm1,pm2,pcm,P2)
% Compute velocities
vMh = jacobian(pMh,q)*dq;
vMt = jacobian(pMt,q)*dq;
vm1 = jacobian(pm1,q)*dq;
vm2 = jacobian(pm2,q)*dq;
vcm = jacobian(pcm,q)*dq;
write_symbolic_term_to_mfile(q,dq,params,vMh,vMt,vm1,vm2,vcm)
% Kinetic energy
K_Mh = (1/2) * Mh * (vMh.'* vMh);
K_Mt = (1/2) * Mt * (vMt.'* vMt); 
K_m1 = (1/2) * m * (vm1.'* vm1); 
K_m2 = (1/2) * m * (vm2.'* vm2); 
K = K_m1 + K_Mh + K_Mt + K_m2;
% Potential energy
V_Mh = pMh(2)*Mh*g;
V_Mt = pMt(2)*Mt*g;
V_m1 = pm1(2)*m*g;
V_m2 = pm2(2)*m*g;
V = V_m1 + V_Mh + V_Mt + V_m2;
% Inertia matrix
D = jacobian(jacobian(K,dq).',dq);
D= simplify(D);
% Coriolis matrix
N = max(size(q));
 syms C
for k = 1:N
    for j = 1:N
        C(k,j) = 0;
        for i = 1:N
            term1 = diff(D(k,j), q(i));
            term2 = diff(D(k,i), q(j));
            term3 = diff(D(i,j), q(k));
            C(k,j) = C(k,j) + 0.5 * (term1 + term2 - term3) * dq(i);
        end
    end
end
C = simplify(C);
% Gravity matrix
G = jacobian(V,q).';
G = simplify(G);
% Control input matrix
B = sym([0 0; 1,0; 0,1]);
% Write matrices to file
write_symbolic_term_to_mfile(q,dq,params,D,C,G,B)
% Impact map
p_e = [p_h; p_v];
qe = [q; p_h; p_v];
dqe = [dq; dp_h; dp_v];
pMh_e = pMh + p_e;
pMt_e = pMt + p_e;
pm1_e = pm1 + p_e;
pm2_e = pm2 + p_e;
P2e = P2 + p_e;
vMh_e = jacobian(pMh_e,qe)*dqe;
vMt_e = jacobian(pMt_e,qe)*dqe;
vm1_e = jacobian(pm1_e,qe)*dqe;
vm2_e = jacobian(pm2_e,qe)*dqe;
K_Mhe = (1/2) * Mh * (vMh_e.'* vMh_e);
K_Mte = (1/2) * Mt * (vMt_e.'* vMt_e);
K_m1e = (1/2) * m * (vm1_e.'* vm1_e);
K_m2e = (1/2) * m * (vm2_e.'* vm2_e);
Ke = K_m1e + K_Mhe + K_Mte + K_m2e;


De = jacobian(jacobian(Ke,dqe).',dqe); De = simplify(De);
E = jacobian(P2e,qe); E= simplify(E);
dY_dq = jacobian(pMh_e,q);
% Write impact map matrices to file
write_symbolic_term_to_mfile(q,dq,params,De,E,dY_dq)
% Controller
fx = [dq; D\(-C*dq - G)];
gx = [zeros(3, 2); D\B];  % Use gx instead of B_ss
syms s delq
%s = (q1 - q1_plus)/delq; 
%delq = q1_minus - q1_plus 
%ds/dt = dq1/delq; ds/dq1 = 1/delq;

syms a21 a22 a23 a24 a25 
syms a31 a32 a33 a34 a35

a2 = [a21 a22 a23 a24 a25];
a3 = [a31 a32 a33 a34 a35];
M = 4;

b2 = 0; b3 = 0;
for k = 0:M
    b2 = b2 + a2(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

for k = 0:M
    b3 = b3 + a3(1,k+1)*(factorial(M)/(factorial(k)*factorial(M-k)))*s^k*(1-s)^(M-k);
end

% Defining outputs

h = [q2 - b2; q3 - b3];

% y_dot = Lfh = dh/dx*fx - independent of gx*u since relative degree is 2
% However, h is a function of (s,q2,q3), not q1 directly, so the following
% is used:
% dh/dq1 = dh/ds*ds/dq1 = dh/ds*1/delq
%
% Temporary variable that multiples the 1st column with 1/delq
temp = sym(eye(6)); temp(1)  = 1/delq;

Lfh = jacobian(h,[s; q2; q3; dq])*temp*fx;

dLfh = jacobian(Lfh,[s; q2; q3; dq])*temp;

% Write matrix used in feedback linearization - d/dx(Lfh) to file
% Inputs:
%       s = (q1 - q1_plus)/delq: gait timing variable
%       delq = q1_minus - q1_plus: difference in cyclic variable during gait 
%       dq1
%       params: 
%       a2: bezier coefficents (1st - 5th) for q2
%       a3: bezier coefficents (1st - 5th) for q3
%
% Outputs:
%       dLfh: partial of Lfh, to be used to compute L2fh and LgLfh
%
write_symbolic_term_to_mfile([s,delq],dq1,[a2,a3],dLfh);



%-------------------------------------------------------------------------%
%%%% For Zero Dynamics

db_ds2 = 0;
for k = 0:M-1
    db_ds2 = db_ds2 + (a2(1,k+2)-a2(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

db_ds3 = 0;
for k = 0:M-1
    db_ds3 = db_ds3 + (a3(1,k+2)-a3(1,k+1))*(factorial(M)/(factorial(k)*factorial(M-k-1)))*s^k*(1-s)^(M-k-1);
end

partial_db_ds2 = jacobian(db_ds2,s)*dq1/delq;

partial_db_ds3 = jacobian(db_ds3,s)*dq1/delq;

beta1 = [partial_db_ds2; partial_db_ds3]*dq1/delq;

beta2 = jacobian(K,dq1);

write_symbolic_term_to_mfile(s,[dq1, delq],[a2, a3],beta1)

write_symbolic_term_to_mfile(q,dq,params,beta2)
%}