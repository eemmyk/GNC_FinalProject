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

nu2_i = mod(nuSolutions_i + 2*pi, 2*pi)

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
