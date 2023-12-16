%% Shape-Based Approach to Low-Thrust Rendezvous Trajectory Design

clear;
close all;

t = linspace(0, 100, 1000);

N = 2;
% phi = 0;
% 
theta_0 = 0;
theta_tilde = 0;

%The total transfer angle is represented by:
theta_f = 2*pi * N + theta_tilde;
%theta = linspace(theta_0, theta_f, 1000);

mju_Earth = 3.986004418*10^14;
mju_Sun = 1.32712440018*10^20;

mju = mju_Earth;

%initial and target orbit parameters
r1 = 6800*1000;%100 * 10^9;
e1 = 0.02;
%It's actually nu not E
E1 = theta_0;
gamma1 = asin(e1 * sin(E1) / (1-e1*cos(E1)));
a1 = r1 / (1 - e1 * cos(E1));
%theta2_dot = d/dt acos((1 - r/a)/e); -->
v1 = sqrt(2*mju / r1 - mju/a1)
theta1_dot = v1/(a1*sqrt(-(r1^2/a1^2) + 2*r1/a1 + e1^2 - 1))

r2 = 15000*1000;%150 * 10^9;
e2 = 0.01;
E2 = theta_f;
gamma2 = asin(e2 * sin(E2) / (1-e2*cos(E2)));
a2 = r2 / (1 - e2* cos(E2));
%theta2_dot = d/dt acos((1 - r/a)/e); -->
v2 = sqrt(2*mju / r2 - mju/a2)
theta2_dot = v2/(a2*sqrt(-(r2^2/a2^2) + 2*r2/a2 + e2^2 -1))

%%
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

syms a b c d f theta;


theta1_dot = sqrt(mju/r1^4) / ((1/r1) + 2*c);
theta2_dot = sqrt(mju/r2^4) / ((1/r2) + 2*c + 6*d*theta_f + 12*e*theta_f^2 + 20 * f*theta_f^3);

%gamma1 = atan(-r1 * b);
%gamma2 = atan(-r2 * (b + 2*c*theta_f + 3*d*theta_f^2 + 4*e*theta_f^3 + 5*f*theta_f^4));

a = 1/r1;
b = -tan(gamma1) / r1;
c = 1/(2*r1) * (mju / (r1^3 * theta1_dot) - 1);

efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
            -48*theta_f     18*theta_f^2 -2*theta_f^3; 
             20            -8*theta_f     theta_f^2];

efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
            -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
            mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];

efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;

e = efg(1);
f = efg(2);
g = efg(3);

%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);

r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
gamma = atan(-r * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4));

f_Theta_dot = sqrt((mju./r.^4)) ./ (1./r + 2*c + 6*d*theta + 2*e*theta.^2 + 20*f*theta.^3 + 30*g*theta.^4);
f_T_a = -mju ./ (2 * r.^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta.^2 + 120*g*theta.^3 - tan(gamma)./r) ./ (1./r + 2*c + 6*d*theta + 2*e*theta.^2 + 20*f*theta.^3 + 30*g*theta.^4).^2;

megaFunction = f_T_a / f_Theta_dot;

costFunction = @(theta) megaFunction;

deltaV = integral(costFunction, theta_0, theta_f);

x = cos(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
y = sin(theta) ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);
r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5);

figure;
hold on;

plot(x, y);

%plot(cos(theta))

%plot(r)
