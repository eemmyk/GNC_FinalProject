% Constants Section

earthRadius = 6.371e6; % m
orbitalAltitude_inner = 2e6; % m
orbitalAltitude_outer = 42.164e6; % m

r1 = orbitalAltitude_inner + earthRadius;
r2 = orbitalAltitude_outer + earthRadius;

% The radius of the orbit can be represented as:
% r(theta) 1/(a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5)

% N is the number of rotations before the orbit settles. 

N = 4;
theta_0 = 0;
theta_tilde = pi/2;

%The total transfer angle is represented by:
theta_f = 2*pi * N + theta_tilde;
theta = linspace(theta_0, theta_f, 1e4);

a = 1/r1;
b = -1e-9;
c = -1e-10;
d = 1e-11;
e = 1e-12;
f = 1e-13;

r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);

% gamma is the flight path angle. The flight path angle is the angle
% between the tangent to Earth and the velocity vector of the spacecraft

% r_dot = -r.^2 * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4)*theta_dot;

% Although this equation is defined, it cannot be used because theta_dot is
% at this point unknown

% tan(gamma1) = -r1 * b
% tan(gamma2) = -r2 * (b + 2*c*theta_f + 3*d*theta_f^2 + 4*e*theta_f^3 + 5*f*theta_f^4)






%theta = theta_0 = 0 yields:
%r1 = 1/a

%theta = theta_f yields:
%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);

% The time derivative of the radius equation yields
% 
% From this we can get the flight






%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);


x = r.*cos(theta);
y = r.*sin(theta);


figure;
hold on;

plot(x, y);