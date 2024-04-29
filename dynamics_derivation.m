D = [Mt*l^2 + Mh*r^2 + Mt*r^2 + (5*m*r^2)/4 + (m*r^4)/4 - m*r^2*cos(q2) - 2*Mt*l*r*cos(2*q1 + q3), -(m*r^2*(2*cos(q2) - 1))/4, Mt*l*(l - r*cos(2*q1 + q3)); -(m*r^2*(2*cos(q2) - 1))/4, (m*r^2)/4, 0; Mt*l*(l - r*cos(2*q1 + q3)), 0, Mt*l^2];
C = [(dq2*m*r^2*sin(q2))/2 + 2*Mt*dq1*l*r*sin(2*q1 + q3) + Mt*dq3*l*r*sin(2*q1 + q3), (m*r^2*sin(q2)*(dq1 + dq2))/2, Mt*l*r*sin(2*q1 + q3)*(dq1 + dq3); -(dq1*m*r^2*sin(q2))/2, 0, 0; Mt*dq1*l*r*sin(2*q1 + q3), 0, 0];
G = [g*m*((r*sin(q1 + q2))/2 - r*sin(q1)) - Mt*g*(l*sin(q1 + q3) + r*sin(q1)) - Mh*g*r*sin(q1) - (g*m*r^2*sin(q1))/2; (g*m*r*sin(q1 + q2))/2; -Mt*g*l*sin(q1 + q3)];
B = [0, 0; 1, 0; 0, 11];
