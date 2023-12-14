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
% From this we can get the flight




a = 1/r1;
b = -1e-9;
c = -1e-10;
d = 1e-11;
e = 1e-12;
f = 1e-13;

%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);


x = cos(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
y = sin(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);

figure;
hold on;

plot(x, y);

%plot(cos(theta))

%plot(r)
