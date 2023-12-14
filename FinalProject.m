t = linspace(0, 100, 1000);

N = 2;
% phi = 0;
% 
theta_0 = 0;
theta_tilde = 0;

%The total transfer angle is represented by:
theta_f = 2*pi * N + theta_tilde;
theta = linspace(theta_0, theta_f, 1000);

r1 = 100 * 10^6;
r2 = 150 * 10^6;

% The radius of the orbit can be represented as:
% r(theta) 1/(a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5)

%theta = theta_0 = 0 yields:
%r1 = 1/a

%theta = theta_f yields:
%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);

% The time derivative of the radius equation yields
% r_dot = -r.^2 * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4)*theta_dot
% From this we can get the flight-path angle gamma may be found
% tan(gamma) = r_dot / (r * theta_dot) = -r * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4)


%The first and final flight path angles gamma1 and gamma2 can be solved from this equation:
%tan(gamma1) = -r1 * b
%tan(gamma2) = -r2 * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4)

%Additionally the polynomial must satisfy the equations of motion
%r_dot_dot = r*theta_dot^2 + mju/r^2 = Ta * sin(alpha)
%1 / r d/dt (r^2*theta_dot) = Ta * cos(alpha)

%Ta is the thrust acceleration and alpha is the thrust angle

%r*theta_dot_dot + 2*r*tan(gamma) = Ta * cos(alpha)

%A LOT OF STUFF HERE:

syms a b c d e f g;


a = 1/r1;
b = -tan(gamma1) / r1;
c = 1/(2*r1) * (mju / (r1^3 * theta_dot_1) - 1);




d = 1e-11;
e = 1e-12;
f = 1e-13;

%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);


f_Theta_a = sqrt((mju/r^4)) / (1/r + 2*c + 6*d*theta + 2*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)
f_T_a = -mju / (2 * r^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta^2 + 120*g*theta^3 - tan(gamma)/r) / (1/r + 2*c + 6*d*theta + 2*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)^2 

costFunction = @(theta) f_T_a / f_Theta_dot;

deltaV = integral(costFunction, theta_0, theta_f);

x = cos(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
y = sin(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);

figure;
hold on;

plot(x, y);

%plot(cos(theta))

%plot(r)
