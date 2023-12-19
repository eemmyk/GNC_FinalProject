clc
clear

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
theta = linspace(theta_0, theta_f, 1e2);

transferTime = 23*60*60; % seconds

%..........................................................................
% The "I don't get how these are supposed to be defined" section
%..........................................................................

theta_dot = theta_f/orbitalPeriod; % rads/s
gamma1 = pi/20;
gamma2 = -pi/20;
d = 0;

%..........................................................................
% Coefficients Definitions
%..........................................................................

% The constants a, b, and c are defined explicitly

a = 1/r1;
b = -tan(gamma1)/r1;
c = (1/2*r1)*((mu/((r1^3)*theta_dot^2)) - 1);

% Next, e, f, and g are found through use of the syms package

syms theta_f_ theta_dot_ r2_ a_ b_ c_ d_ e_ f_ g_ gamma2_ mu_ r_ theta_ tf_

efgMatrix1 = [30*theta_f_^2, -10*theta_f_^3, theta_f_^4;
              -48*theta_f_, 18*theta_f_^2, -2*theta_f_^3;
              20, -8*theta_f_, theta_f_^2];
efgMatrix2 = [(1/r2_) - (a_ + b_*theta_f_ + c_*theta_f_^2 + d_*theta_f_^3);
              -(tan(gamma2_)/r2_) - (b_ + 2*c_*theta_f_ + 3*d_*theta_f_^2);
              (mu_/((r2_^4)*theta_dot_^2)) - ((1/r2_) + 2*c_ + 6*d_*theta_f_)];
efgSolution_sym = (1/(theta_f_^6))*efgMatrix1*efgMatrix2;

% efgSolution = subs(efgSolution_sym, [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_], ...
%                    [theta_f, theta_dot, r2, a, b, c, d, gamma2, mu]);
% 
% efgSolution = double(efgSolution);
% 
% e = efgSolution(1);
% f = efgSolution(2);
% g = efgSolution(3);

%..........................................................................
% Finding the transfer time calculated by the constants
%..........................................................................

% This section is deprecated because it was attempting analytical
% integration

% rFunc = 1 ./ (a_ + b_*theta_ + c_*theta_.^2 + d_*theta_.^3 + e_*theta_.^4 ...
%               + f_*theta_.^5 + g_*theta_.^6);
% 
% tfFunction = sqrt(((rFunc^4)/mu_)*((1/rFunc) + 2*c_ + 6*d_*theta_ + 12*e_*theta_^2 ...
%                   + 20*f_*theta_^3 + 30*g_*theta_^4));
% 
% t_f_ = int(tfFunction, theta_, [0, tf_]);

% rFunc = 1 ./ (a + b*thetaInt + c*thetaInt.^2 + d*thetaInt.^3 + e*thetaInt.^4 ...
%               + f*thetaInt.^5 + g*thetaInt.^6);

% debugTest = real(sqrt(((rFunc.^4)/mu).*((1./rFunc) + 2*c + 6*d*theta + 12*e*theta.^2 ...
%                    + 20*f*theta.^3 + 30*g*theta.^4)));

% integrateTransferTime = @(thetaInt) real(sqrt((((1 ./ (a + b*thetaInt + c*thetaInt.^2 + d*thetaInt.^3 + e*thetaInt.^4 ...
%               + f*thetaInt.^5 + g*thetaInt.^6)).^4)/mu).*((1./(1 ./ (a + b*thetaInt + c*thetaInt.^2 + d*thetaInt.^3 + e*thetaInt.^4 ...
%               + f*thetaInt.^5 + g*thetaInt.^6))) + 2*c + 6*d*thetaInt + 12*e*thetaInt.^2 ...
%                   + 20*f*thetaInt.^3 + 30*g*thetaInt.^4)));
% 
% tf_integralValue = integral(integrateTransferTime, 0, theta_f);
% 
% dFinder = fzero(integrateTransferTime - transferTime, 0);

bestGuess = 1e10;

firstBound = 1e20;
boundStep = 1e20;
secondBound = 1e25;

% dValue = fzero(@dOptimiser, 0);
dValue = fzero(@dOptimiser, 1e-15);

% 
% dOptimiser(1e20)

function dOptimised = dOptimiser(dGuess)
    earthRadius = 6.371e6; % m
    orbitalAltitude_inner = 2e6; % m
    orbitalAltitude_outer = 42.164e6; % m
    G = 6.6743e-11; % gravitational parameter Nm^2/kg^2
    M_earth = 5.972e+24; % mass of Earth, kg
    mu = G*M_earth;
    
    r1 = orbitalAltitude_inner + earthRadius;
    r2 = orbitalAltitude_outer + earthRadius;
    N = 1; % number of loops
    theta_0 = 0; % r1 angle
    theta_tilde = pi/2; % angle between r1 and r2
    
    semiMajor = (r1+r2)/2;
    orbitalPeriod = 2*pi*((semiMajor^3)/mu);

    theta_f = 2*pi * N + theta_tilde;
    theta_dot = theta_f/orbitalPeriod; % rads/s
    gamma1 = pi/20;
    gamma2 = -pi/20;
    a = 1/r1;
    b = -tan(gamma1)/r1;
    c = (1/2*r1)*((mu/((r1^3)*theta_dot^2)) - 1);
    % Definitions since MATLAB hates global variables
    syms theta_f_ theta_dot_ r2_ a_ b_ c_ d_ e_ f_ g_ gamma2_ mu_ r_ theta_ tf_
    symbolicInputs = [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_];
    
    efgMatrix1 = [30*theta_f_^2, -10*theta_f_^3, theta_f_^4;
                  -48*theta_f_, 18*theta_f_^2, -2*theta_f_^3;
                  20, -8*theta_f_, theta_f_^2];
    efgMatrix2 = [(1/r2_) - (a_ + b_*theta_f_ + c_*theta_f_^2 + d_*theta_f_^3);
                  -(tan(gamma2_)/r2_) - (b_ + 2*c_*theta_f_ + 3*d_*theta_f_^2);
                  (mu_/((r2_^4)*theta_dot_^2)) - ((1/r2_) + 2*c_ + 6*d_*theta_f_)];
    efgSolution_sym = (1/(theta_f_^6))*efgMatrix1*efgMatrix2;
    transferTime = 23*60*60; % seconds
    efgSolution = subs(efgSolution_sym, symbolicInputs, [theta_f, theta_dot, r2, a, b, c, dGuess, gamma2, mu]);
    efgSolution = double(efgSolution);
    e = efgSolution(1);
    f = efgSolution(2);
    g = efgSolution(3);
    integrateTransferTime = @(thetaInt) abs(sqrt((((1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
              + f*thetaInt.^5 + g*thetaInt.^6)).^4)/mu).*((1./(1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
              + f*thetaInt.^5 + g*thetaInt.^6))) + 2*c + 6*dGuess*thetaInt + 12*e*thetaInt.^2 ...
                  + 20*f*thetaInt.^3 + 30*g*thetaInt.^4)));
    tf_integralValue = integral(integrateTransferTime, 0, theta_f);
    dOptimised = abs(tf_integralValue - transferTime);
end


% for dGuess = firstBound:boundStep:secondBound
%     efgSolution = subs(efgSolution_sym, [theta_f_, theta_dot_, r2_, a_, b_, c_, d_, gamma2_, mu_], ...
%                    [theta_f, theta_dot, r2, a, b, c, dGuess, gamma2, mu]);
%     efgSolution = double(efgSolution);
%     e = efgSolution(1);
%     f = efgSolution(2);
%     g = efgSolution(3);
%     integrateTransferTime = @(thetaInt) abs(sqrt((((1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
%               + f*thetaInt.^5 + g*thetaInt.^6)).^4)/mu).*((1./(1 ./ (a + b*thetaInt + c*thetaInt.^2 + dGuess*thetaInt.^3 + e*thetaInt.^4 ...
%               + f*thetaInt.^5 + g*thetaInt.^6))) + 2*c + 6*dGuess*thetaInt + 12*e*thetaInt.^2 ...
%                   + 20*f*thetaInt.^3 + 30*g*thetaInt.^4)));
%     tf_integralValue = integral(integrateTransferTime, 0, theta_f);
%     minimisationValue = abs(tf_integralValue - transferTime);
%     if minimisationValue <= bestGuess
%         bestGuess = minimisationValue;
%         bestdGuess = dGuess;
%     end
%     fprintf('Running: %2f completed.\n', (((dGuess - firstBound)/(secondBound-firstBound)))*100)
% end




%..........................................................................
% Plotting Section
%..........................................................................

% r = 1 ./ (a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5 + g*theta.^6);
% 
% x = r.*cos(theta);
% y = r.*sin(theta);
% 
% plot(x, y)
% axis equal

