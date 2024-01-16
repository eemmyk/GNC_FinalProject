clc
clear

global theta_f_ theta_dot_ r2_ a_ b_ c_ d_ e_ f_ g_ gamma2_ mu_ r_ theta_ tf_
global earthRadius orbitalAltitude_inner orbitalAltitude_outer G M_earth mu
global semiMajor_inner orbitalPeriod_inner r1 e1 theta_dot1
global semiMajor_outer orbitalPeriod_outer r2 e2 theta_dot2
global N theta_0 theta_tilde theta_f theta transferTime
global a b c e f g gamma1 gamma2
global efgMatrix1 efgMatrix2 efgSolution_sym

%..........................................................................
% Constants Section
%..........................................................................

earthRadius = 6.371e6; % m
orbitalAltitude_inner = 2e6; % m
orbitalAltitude_outer = 42.164e6; % m
G = 6.6743e-11; % gravitational parameter Nm^2/kg^2
M_earth = 5.972e+24; % mass of Earth, kg
mu = G*M_earth;

%..........................................................................
% Orbital Parameters
%..........................................................................

% Transfer Path

N = 1; % number of transfer loops
theta_0 = 0; % r1 angle
theta_tilde = pi/2; % angle between r1 and r2
theta_f = 2*pi * N + theta_tilde;
theta = linspace(theta_0, theta_f, 1e4);

transferTime = 23*60*60; % seconds

% Inner Orbit

semiMajor_inner = orbitalAltitude_inner + earthRadius;
orbitalPeriod_inner = 2*pi*sqrt((semiMajor_inner^3)/mu);
e1 = 0.4;
r1_orbit = (semiMajor_inner*(1 - e1^2))./(1 + e1*cos(theta));
r1 = r1_orbit(1);
theta_dot1 = 1/sqrt((semiMajor_inner^3)/mu);

% Outer Orbit

semiMajor_outer = orbitalAltitude_outer + earthRadius;
orbitalPeriod_outer = 2*pi*sqrt((semiMajor_outer^3)/mu);
e2 = 0.1;
r2_orbit = (semiMajor_outer.*(1 - e2.^2))./(1 + e2.*cos(theta));
r2 = semiMajor_outer;
theta_dot2 = 1/sqrt((semiMajor_outer^3)/mu);

%..........................................................................
% The "I don't get how these are supposed to be defined" section
%..........................................................................

gamma1 = pi/20;
gamma2 = -pi/20;

%..........................................................................
% Coefficients Definitions
%..........................................................................

% The constants a, b, and c are defined explicitly

a = 1/r1;
b = -tan(gamma1)/r1;
c = (1/(2*r1))*((mu/((r1^3)*theta_dot1^2)) - 1);

% Next, e, f, and g are found through use of the syms package

% syms theta_f_ theta_dot_ r2_ a_ b_ c_ d_ e_ f_ g_ gamma2_ mu_ r_ theta_ tf_
% 
% efgMatrix1 = [30*theta_f_^2, -10*theta_f_^3, theta_f_^4;
%               -48*theta_f_, 18*theta_f_^2, -2*theta_f_^3;
%               20, -8*theta_f_, theta_f_^2];
% efgMatrix2 = [(1/r2_) - (a_ + b_*theta_f_ + c_*theta_f_^2 + d_*theta_f_^3);
%               -(tan(gamma2_)/r2_) - (b_ + 2*c_*theta_f_ + 3*d_*theta_f_^2);
%               (mu_/((r2_^4)*theta_dot_^2)) - ((1/r2_) + 2*c_ + 6*d_*theta_f_)];
% efgSolution_sym = (1/(2*theta_f_^6))*efgMatrix1*efgMatrix2;
% 
% efgSolution = subs(efgSolution_sym, [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_], ...
%                    [theta_f, theta_dot2, r2, a, b, c, d, gamma2, mu]);
% 
% efgSolution = double(efgSolution);
% 
% e = efgSolution(1);
% f = efgSolution(2);
% g = efgSolution(3);

%..........................................................................
% Finding the transfer time calculated by the constants
%..........................................................................

d = fminsearch(@dOptimiser, 1e-12);

%..........................................................................
% Plotting Section
%..........................................................................

r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5 + g*theta.^6);

x = r.*cos(theta);
y = r.*sin(theta);

hold on
plot(x, y, 'LineStyle', '--', 'DisplayName', 'Trajectory', 'Color', 'k')
plot(0, 0, 'Marker', 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'Earth', 'MarkerEdgeColor', 'b', 'LineStyle', 'none')

plot(cos(theta).*r1_orbit, sin(theta).*r1_orbit, 'Color', 'k', 'DisplayName', 'Inner Orbit')

plot(cos(theta).*r2_orbit, sin(theta).*r2_orbit, 'Color', 'k', 'DisplayName', 'Outer Orbit')

hold off
xlabel('Distance From Earth (m)')
ylabel('Distance From Earth (m)')
legend
axis equal

%..........................................................................
% Functions Below
%..........................................................................

function dOptimised = dOptimiser(dGuess)
    global theta_f_ theta_dot_ r2_ a_ b_ c_ d_ e_ f_ g_ gamma2_ mu_ r_ theta_ tf_
    global earthRadius orbitalAltitude_inner orbitalAltitude_outer G M_earth mu
    global semiMajor_inner orbitalPeriod_inner r1 e1 theta_dot1
    global semiMajor_outer orbitalPeriod_outer r2 e2 theta_dot2
    global N theta_0 theta_tilde theta_f theta transferTime
    global a b c e f g gamma1 gamma2
    global efgMatrix1 efgMatrix2 efgSolution_sym

    efgSolution = subs(efgSolution_sym, [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_], [theta_f, theta_dot2, r2, a, b, c, dGuess, gamma2, mu]);
    efgSolution = double(efgSolution);
    e = efgSolution(1);
    f = efgSolution(2);
    g = efgSolution(3);
    integrateTransferTime = @(thetaInt) abs(sqrt((((1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
              + f*thetaInt.^5 + g*thetaInt.^6)).^4)/mu).*((1./(1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
              + f*thetaInt.^5 + g*thetaInt.^6))) + 2*c + 6*dGuess*thetaInt + 12*e*thetaInt.^2 ...
                  + 20*f*thetaInt.^3 + 30*g*thetaInt.^4)));
    tf_integralValue = integral(integrateTransferTime, 0, theta_f);
    dOptimised = abs(transferTime - tf_integralValue);
end






