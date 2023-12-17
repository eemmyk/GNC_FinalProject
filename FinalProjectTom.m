%..........................................................................
% Constants Section
%..........................................................................

earthRadius = 6.371e6; % m
orbitalAltitude_inner = 2e6; % m
orbitalAltitude_outer = 42.164e6; % m
G = 6.6743e-11; % gravitational parameter Nm^2/kg^2
M_earth = 5.972e+24; % mass of Earth, kg
mu = G*M_earth;

r1 = orbitalAltitude_inner + earthRadius;
r2 = orbitalAltitude_outer + earthRadius;

%..........................................................................
% Orbital Parameters
%..........................................................................

N = 1; % number of loops
theta_0 = 0; % r1 angle
theta_tilde = pi/2; % angle between r1 and r2

semiMajor = (r1+r2)/2;
orbitalPeriod = 2*pi*((semiMajor^3)/mu);

theta_f = 2*pi * N + theta_tilde;
theta = linspace(theta_0, theta_f, 1e6);

transferTime = 23*60*60; % seconds

%..........................................................................
% The "I don't get how these are supposed to be defined" section
%..........................................................................

theta_dot = (2*pi)/orbitalPeriod; % rads/s
gamma1 = pi/20;
gamma2 = -pi/20;
d = 1e+10;

%..........................................................................
% Coefficients Definitions
%..........................................................................

% The constants a, b, and c are defined explicitly

a = 1/r1;
b = -tan(gamma1)/r1;
c = (1/2*r1)*((mu/((r1^3)*theta_dot^2)) - 1);

% Next, e, f, and g are found through use of the syms package

syms theta_f_ theta_dot_ r2_ a_ b_ c_ d_ gamma2_ mu_

efgMatrix1 = [30*theta_f_^2, -10*theta_f_^3, theta_f_^4;
              -48*theta_f_, 18*theta_f_^2, -2*theta_f_^3;
              20, -8*theta_f_, theta_f_^2];
efgMatrix2 = [(1/r2_) - (a_ + b_*theta_f_ + c_*theta_f_^2 + d_*theta_f_^3);
              -(tan(gamma2_)/r2_) - (b_ + 2*c_*theta_f_ + 3*d_*theta_f_^2);
              (mu_/((r2_^4)*theta_dot_^2)) - ((1/r2_) + 2*c_ + 6*d_*theta_f_)];
efgSolution_sym = (1/(theta_f_^6))*efgMatrix1*efgMatrix2;

efgSolution_sym = subs(efgSolution_sym, [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_], ...
                   [theta_f, theta_dot, r2, a, b, c, d, gamma2, mu]);

efgSolution = double(efgSolution_sym);

e = efgSolution(1);
f = efgSolution(2);
g = efgSolution(3);

%..........................................................................
% Plotting Section
%..........................................................................

r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5 + g*theta.^6);

x = r.*cos(theta);
y = r.*sin(theta);

plot(x, y)
axis equal

% Note to self:
% Where I left off, I was about to put in the integral which autocalculates
% the constant 'd'. I got confused because the integral requires all the
% other constants, but e f and g are dependent on d. This is a tad
% complicated but I think the point is to specify a loop where we calculate
% for each loop of d and after each integral, check if the difference
% between the two values is smallest. Also, it's possible that theta_dot
% might be easily calculatable from the specified transfer time.

