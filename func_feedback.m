% Compute control action using feedback linearization
%
% Inputs:
%       x: states
%           q1
%           q2
%           q3
%           dq1
%           dq2
%           dq3
%       alpha: Bezier coefficients for q2 and q3
%           alpha2 (1st to 5th coefficients)
%           alpha3 (1st to 5th coefficients)
%       s_params: for gait timing
%           q1_min
%           delq
%
% Outputs:
%       u: control action
%

function u = func_feedback(x,alpha,s_params)
% gains
kp1 = 700;
kp2 = 700;
kd1 = 25;
kd2 = 25;

% Seperating inputs
q = x(1:3);
dq = x(4:6);

% Get model parameters
[r,m,Mh,Mt,l,g] = func_model_params;
params = [r,m,Mh,Mt,l,g];

% Seperate Bezier coefficients
alpha2 = alpha(1:5);
alpha3 = alpha(6:10);

% Seperate s_params
q1_min = s_params(1);
q1_max = s_params(2);
delq = q1_max - q1_min;

% Gait timing variable
% Inputs:
%       q1
%       q1_min
%       delq
s = func_gait_timing(q(1), q1_min, q1_max);

% Building output function
%       y = h(x) = Hq - hd

% Get D,C,G,B matrices
% Inputs:
%       q = [q1, q2, q3]
%       dq = [dq1, dq2, dq3]
%       params = [r,m,Mh,Mt,l,g]
[D,C,G,B] = func_compute_D_C_G_B(q,dq,params);

% Defining fx and gx
fx = [dq; inv(D) * (-C*dq-G)]; % 6x1
gx = [zeros(3,2); inv(D)*B]; % 6x2

% find y
b2 = bezier(s, 4, alpha2);
b3 =  bezier(s, 4, alpha3);
h = [q(2) - b2; q(3) - b3];
y = h; % 2x1

dh_dx = [ 0, 1, 0, 0, 0, 0; 
          0, 0, 1, 0, 0, 0]; % 2x6
dh_dx(1,1) = -d_ds_bezier(s, 4, alpha2)/delq;
dh_dx(2,1) = -d_ds_bezier(s, 4, alpha3)/delq;  
q = [q(2); q(3)];

Lfh = dh_dx*fx;
dy = Lfh; % 2x1x

%%%% PD controller
Kp = [kp1,0; 0,kp2];
Kd = [kd1,0; 0,kd2];
epsilon = 0.5;
v = (1/epsilon^2)*Kp*y + (1/epsilon)*Kd*dy; % 2x1

%%%% Feedback linearization:
dLfh = func_compute_dLfh([s,delq],dq(1),[alpha2,alpha3]);
L2fh =  dLfh*fx; % 2x1
LgLfh = dLfh*gx; % 2x2

% Control action that uses feedback linearization with PD controller
u = -pinv(LgLfh) * (v + L2fh); % 2x1
end