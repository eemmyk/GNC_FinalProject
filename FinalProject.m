%% Shape-Based Approach to Low-Thrust Rendezvous Trajectory Design

clear;
close all;

%Things to change to test stuff out:
% tf, desired time of flight
% N, number of rotations around central body
% theta_tilde, angle between initial and end position
% Any of the orbital parameters
% limits for d
% And stuff

global tf  theta_f  theta_0  integralApproximationSteps;

%% Estimation values and set parameters

%Gravitational parameter
%mju = 3.986004418*10^14; %Earth
mju = 1.32712440018*10^20; %Sun

%First orbit semimajor axis
%a_initial= 6800*1000;%500 * 10^9;
a_initial = 1050 * 10^9;
%Second orbit semimajor axis
%a_final = 42164*1000;
a_final = 2125*10^9;
%Number of rotations around central body
N = 0;

%Solve TOF or optimize deltaV
optimizeDV = 1; % 1 --> optimize DV, 0--> solve exact TOF

% 1 + 2N hohmann transfer orbits in time
TOF_estimation = 0.5*(1+N)*pi*sqrt((a_initial+a_final)^3/(8*mju));

%tf = 86400 * 365.25 * 0.5;
%ratio = tf/TOF_estimation;

tf = TOF_estimation;

integralApproximationSteps = 10;

%%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

%randomized variables if wanting to use them

% Target will be at random point along 2nd orbit
currentTime = 0; %rand() * 2*pi*sqrt(a_final ^3 / mju);

% Initial burn will be performed at a random point along the first orbit
nu_0 = 0; %rand() * 2 * pi;

%% Calculating the orbital parameters

%initial orbit parameters
a1 = a_initial;
e1 = 0.4;
%Argument of perigee
omega1 = pi/4;

%initial maneuver angle
nu1 = omega1 + nu_0;
%In reference coords
theta1 = nu1;
gamma1 = asin(e1 * sin(nu_0) / sqrt(1+2*e1*cos(nu_0) + e1^2));
p1 = a1 * (1-e1^2);
r1 = p1 / (1+e1*cos(nu_0));
theta1_dot = sqrt(mju/a1^3) * a1^2/r1^2 * sqrt(1-e1^2);

%Target orbit parameters
%Some are known, while others are calculated from desired tof and initial
%conditions
a2 = a_final;
e2 = 0.7;
%Argument of perigee
omega2 = -pi/2;
%Time of last perigee pass for object 2
Tp2 = 0;

%Trying this one instead
n2 = sqrt(mju/a2^3);
P2 = 2*pi/n2;
T = mod(currentTime, P2) - Tp2;

% nT = E - e*sin(E)
% E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));

syms nu_time

E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
nuSolver = n2*T==pi+E-e2*sin(E);
nuSolutions_i = vpasolve(nuSolver, nu_time);
nuSolutions_i = double(nuSolutions_i);

nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);

%Now we can continue by calculating the true anomaly when the spacecraft
%reaches orbit 2. At currentTime + tf

finalTime = currentTime + tf;
T = mod(finalTime, P2) - Tp2;

% nT = E - e*sin(E)
% E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));

syms nu_time

E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
nuSolver = n2*T==pi+E-e2*sin(E);
nuSolutions_f = vpasolve(nuSolver, nu_time);
nuSolutions_f = double(nuSolutions_f);

nu2 = mod(nuSolutions_f + 2*pi, 2*pi);
%In reference coords
theta2 = nu2 + omega2;
theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);

gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
p2 = a2 * (1-e2^2);
r2 = p2 / (1+e2*cos(nu2));
r2_i = p2 / (1+e2*cos(nu2_i));
theta2_dot = sqrt(mju/a2^3)*a2^2/r2^2 * sqrt(1-e2^2);


%The total transfer angle is represented by:

theta_f = 2*pi * N + theta_tilde;

%% Solving the coefficients
syms d theta;

a = 1/r1;
b = -tan(gamma1) / r1;
c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);

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

r = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);
gamma = atan(-r * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4 + 6*g*theta^5));

global timeFunction thrustFunction thetaDotFunction;


thetaDotFunction = sqrt((mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
thrustFunction = -mju / (2 * r^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta^2 + 120*g*theta^3 - tan(gamma)/r) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)^2;

%Initial guess for d coefficient:
d_guess = 1e-12;

%% TOF optimization and result plotting
if optimizeDV == 0
    %theta_cross = fminsearch(@thetaSquareFunc, [0, 0]);
    
    timeFunction = sqrt((r^4/mju) * (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
    
    opt = optimset('TolFun', 1e-6);

    %TOF solution
    d_optimized_tof = fzero(@transferTimeOptimization, d_guess, opt);
    

    %Plot results of TOF solved trajectory
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a1,a2) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    a_n = a;
    b_n = b;
    c_n = c;
    
    theta_n = double(subs(theta, theta, linspace(theta_0, theta_f, 1000)));
    e_n = double(subs(e, d, d_optimized_tof));
    f_n = double(subs(f, d, d_optimized_tof));
    g_n = double(subs(g, d, d_optimized_tof));
    d_n = double(subs(d, d, d_optimized_tof));
    
    r_es = 1 ./ (a_n + b_n*theta_n + c_n*theta_n.^2 + d_n*theta_n.^3 + e_n*theta_n.^4 + f_n*theta_n.^5 + g_n*theta_n.^6);
    x = cos(theta_n+theta1) .* r_es;
    y = sin(theta_n+theta1) .* r_es;
    
    plot(x, y, "Color", [1 0.1 0.1]);
    
    title("TOF optimized trajectory")
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end

%% Delta V optimization and result plotting
if optimizeDV == 1
    d_fuelOptimal_guess = real(d_guess);
    
    opt = optimset('TolFun',1e1);
    d_fuelOptimal = fminsearch(@deltaVOptimization, d_fuelOptimal_guess, opt);
    
    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_fuelOptimal);
    thrustFunction_n = subs(thrustFunction, d, d_fuelOptimal);
    theta_vec = linspace(theta_0, theta_f, integralApproximationSteps);
    jerk = zeros(2, integralApproximationSteps);
    thrust = zeros(2, integralApproximationSteps);
    for i = 1:integralApproximationSteps
        jerk(:,i) = [theta_vec(i); double(subs(jerkFunction_n, theta, theta_vec(i)))];
        thrust(:,i) = [theta_vec(i); double(subs(thrustFunction_n, theta, theta_vec(i)))];
    end
    
    deltaV = trapz(jerk(1, :), abs(jerk(2,:)));
    
    figure;
    plot(thrust(1, :), thrust(2, :));
    
    
    timeFunction_n = subs(timeFunction, d, d_fuelOptimal);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    dV_TOF = integral(timeFunction_nn, theta_0, theta_f);
    
    title("Thrust curve")
    xlabel("theta");
    ylabel("thurst [N]");
    
    nu = linspace(0, 2*pi, 1000);
    orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
    orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];
    
    % Plot results of dV optimized trajectory
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a1,a2) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    a_n = a;
    b_n = b;
    c_n = c;
    
    theta_n = double(subs(theta, theta, linspace(theta_0, theta_f, 1000)));
    e_n = double(subs(e, d, d_fuelOptimal));
    f_n = double(subs(f, d, d_fuelOptimal));
    g_n = double(subs(g, d, d_fuelOptimal));
    d_n = double(subs(d, d, d_fuelOptimal));
    
    r_es = 1 ./ (a_n + b_n*theta_n + c_n*theta_n.^2 + d_n*theta_n.^3 + e_n*theta_n.^4 + f_n*theta_n.^5 + g_n*theta_n.^6);
    x = cos(theta_n+theta1) .* r_es;
    y = sin(theta_n+theta1) .* r_es;
    
    plot(x, y, "Color", [1 0.1 0.1]);
    
    
    title("deltaV optimized trajectory")
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end