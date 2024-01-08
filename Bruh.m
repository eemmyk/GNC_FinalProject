%Define starting true anomaly
%nu2_i = 0;

%Time of last perigee pass for object 2
Tp2 = 0;
currentTime = 90000;

%This can be used to calculate the current true anomaly with:
% M = nT = E - e*sin(E)
% sin(E) = sin(nu2_tf)*sqrt(1-e2^2) / (1+e2*cos(nu2_tf))
% n2 = sqrt(mju/a2^3);
% T = currentTime - Tp2;
% 
% syms nu2_i
% 
% nuSolver_i = n2*T==asin(sin(nu2_i)*sqrt(1-e2^2)/(1+e2*cos(nu2_i)))-e2*sin(nu2_i)*sqrt(1-e2^2) / (1+e2*cos(nu2_i));
% nuSolutions_i = vpasolve(nuSolver_i, nu2_i);
% nuSolutions_i = double(nuSolutions_i);



%Trying this one instead
n2 = sqrt(mju/a2^3);
P2 = 2*pi/n2;
T = mod(currentTime, P2) - Tp2;

% M = nT = E - e*sin(E)
% E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));

% n2 = sqrt(mju/a2^3);
% T = currentTime - Tp2;


syms nu2_i

E = 2*atan(tan((nu2_i-pi)/2)/sqrt((1+e2)/(1-e2)));
nuSolver_i = n2*T==pi+E-e2*sin(E);
nuSolutions_i = vpasolve(nuSolver_i, nu2_i);
nuSolutions_i = double(nuSolutions_i);

nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);

%%
%Solve nu at tf based on kelper's equations:
%Calculating the true anomaly at t = tf
% (M-M0)/n = t-t0
% n = sqrt(mju / a^3)
% M = E-e*sin(E)
% E = acos((e+cos(nu)) / (1+e*cos(nu)))

% sinh(H) = sin(nu) * sqrt(e2^2-1) / (1 + e*cos(nu))

% Xi = e2*sin(nu2_i)*sqrt(e2^2-1)/(1+e2*cos(nu2_i))-asinh(sin(nu2_i)*sqrt(e2^2-1)/(1 + e2*cos(nu2_i)));
% 
% syms nu2_tf
% 
% nuSolver = sqrt(-a2^3 /mju)*(e2*sin(nu2_tf)*sqrt(e2^2-1)/(1 + e2*cos(nu2_tf))-asinh(sin(nu2_tf)*sqrt(e2^2-1)/(1 + e2*cos(nu2_tf)))-Xi) == tf;
% nuSolutions = solve(nuSolver, nu2_tf);
% nuSolutions_n = double(nuSolutions);
% 
% n2 = sqrt(mju/a2^3);
% M0 = acos((e2+cos(nu2_i))/(1+e2*cos(nu2_i)))-e2*(sin(nu2_i)*sqrt(1-e2^2)/(1+e2*cos(nu2_i)));
% 
% syms nu2_tf
% nuSolver = ((acos((e2+cos(nu2_tf))/(1+e2*cos(nu2_tf)))-e2*(sin(nu2_tf)*sqrt(1-e2^2)/(1+e2*cos(nu2_tf))))-M0)/n2 == tf;
% nuSolutions = vpasolve(nuSolver, nu2_tf);
% nuSolutions_n = double(nuSolutions);
% %Since no retrograde orbits are considered, select the positive nu
% nu2 = nuSolutions_n > 0
% 
% % M = 2*pi*(t)

% The radius of the orbit can be represented as:
% r(theta) 1/(a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5)

%theta = theta_0 = 0 yields:
%r1 = 1/a

%theta = theta_f yields:
%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);

% The time derivative of the radius equation yields
% r_dot = -r.^2 * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)*theta_dot
% From this we can get the flight-path angle gamma may be found
% tan(gamma) = r_dot / (r * theta_dot) = -r * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)


%The first and final flight path angles gamma1 and gamma2 can be solved from this equation:
%tan(gamma1) = -r1 * b
%tan(gamma2) = -r2 * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)

%Additionally the polynomial must satisfy the equations of motion
%r_dot_dot = r*theta_dot^2 + mju/r^2 = Ta * sin(alpha)
%1 / r d/dt (r^2*theta_dot) = Ta * cos(alpha)

%Ta is the thrust acceleration and alpha is the thrust angle

%r*theta_dot_dot + 2*r*tan(gamma) = Ta * cos(alpha)

%A LOT OF STUFF HERE: