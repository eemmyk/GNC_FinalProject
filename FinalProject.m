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

global tf theta_f theta_0 intApprox a_initial a_final;
global nu_0 currentTime mju N;

%% Estimation values and set parameters

%Gravitational parameter
%mju = 3.986004418*10^14; %Earth
mju = 1.32712440018*10^20; %Sun

%First orbit semimajor axis
%a_initial= 6800*1000;%500 * 10^9;
a_initial = 10 * (0.01 + 0.9 * rand()) * 150*10^9;
%Second orbit semimajor axis
%a_final = 42164*1000;
a_final = 10 * (0.01 + 0.9 * rand()) * 150*10^9;
%Number of rotations around central body
N = randi(2) - 1;

%Solve TOF, optimize deltaV or both
optimizeDV = 1;
optimizeTOF = 1;

% 1 + 2N hohmann transfer orbits in time
TOF_estimation = 0.5*(1+N)*pi*sqrt((a_initial+a_final)^3/(8*mju));

%tf = 86400 * 365.25 * 0.5;
%ratio = tf/TOF_estimation;

if optimizeTOF == 1
    tf = TOF_estimation;
else
    tf = 0;
end

%Used for dV optimization
tf_guess = 0.001 * TOF_estimation;

intApprox = 10;
plotAccuracy = 1000;

%%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

%randomized variables if wanting to use them

% Target will be at random point along 2nd orbit
currentTime = rand() * 2*pi*sqrt(a_final ^3 / mju);

% Initial burn will be performed at a random point along the first orbit
nu_0 = rand() * 2 * pi;

%% Calculating the orbital parameters
global theta2 omega1 omega2 e1 e2;

%initial orbit parameters
a1 = a_initial;
e1 = rand()^2;
%Argument of perigee
omega1 = rand()*2*pi;

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
e2 = rand()^2;
%Argument of perigee
omega2 = rand()*2*pi;
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

%Propagating the two orbits
nu = linspace(0, 2*pi, 1000);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

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
global d_solution;

d_solution = 1e-12;

%% TOF optimization and result plotting
global timeResult

%Initialize best values
timeResult = Inf;
if optimizeTOF == 1
    %theta_cross = fminsearch(@thetaSquareFunc, [0, 0]);
    
    timeFunction = sqrt((r^4/mju) * (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
    
    opt = optimset('TolFun', 1e-18, 'TolX', 1e-18);

    %TOF solution
    d_optimized_tof = fzero(@transferTimeOptimization, d_solution, opt);
    
    theta_vec = linspace(theta_0, theta_f, intApprox);
    
    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_optimized_tof);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    
    thrustFunction_n = subs(thrustFunction, d, d_solution);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));

    jerk = [theta_vec; jerkFunction_nn(theta_vec)];    
    thrustCurve_tof = [theta_vec; thrustFunction_nn(theta_vec)];    

    deltaV_tof = trapz(jerk(1, :), abs(jerk(2,:)));
    
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
    
    efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
            -48*theta_f     18*theta_f^2 -2*theta_f^3; 
             20            -8*theta_f     theta_f^2];

    efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
            -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
            mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];

    efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;

    a_n = a;
    b_n = b;
    c_n = c;
    e_n = double(subs(efg(1), d, d_optimized_tof));
    f_n = double(subs(efg(2), d, d_optimized_tof));
    g_n = double(subs(efg(3), d, d_optimized_tof));
    d_n = double(subs(d, d, d_optimized_tof));

    theta_n = double(subs(theta, theta, linspace(theta_0, theta_f, plotAccuracy)));
    
    r_es = 1 ./ (a_n + b_n*theta_n + c_n*theta_n.^2 + d_n*theta_n.^3 + e_n*theta_n.^4 + f_n*theta_n.^5 + g_n*theta_n.^6);
    x = cos(theta_n+theta1) .* r_es;
    y = sin(theta_n+theta1) .* r_es;
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    title(sprintf("TOF optimized trajectory\nTarget time: %s\nAchieved time: %s\nRequired deltaV: %.0f", secToTime(tf), secToTime(tf+timeResult), deltaV_tof));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end

%% Delta V optimization and result plotting
if optimizeDV == 1
    
    global interDeltaResult deltaResult theta2_opt r2_opt tof_optimal;
    global gamma2_opt theta2_dot_opt P2_opt;

    %Initialize best values
    deltaResult = Inf;
    interDeltaResult = Inf;

    opt = optimset('TolFun',1e1);
    %tf_fuelOptimal = fmincon(@optimalDVSolver, tf_guess, [], [], [], [], 0, [], [], opt);

    tf_fuelOptimal = fminsearch(@optimalDVSolver, tf_guess, opt);

    thrustFunction_n = subs(thrustFunction, d, d_solution);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));

    theta_vec = linspace(theta_0, theta_f, intApprox);
    thrustCurve_dV = [theta_vec; thrustFunction_nn(theta_vec)];    

    %timeFunction_n = subs(timeFunction, d, d_fuelOptimal);
    %timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    %dV_TOF = integral(timeFunction_nn, theta_0, theta_f);
    
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a1,a2) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
            -48*theta_f     18*theta_f^2 -2*theta_f^3; 
             20            -8*theta_f     theta_f^2];

    efg_Mat_2 = [1/r2_opt - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
            -tan(gamma2_opt)/r2_opt - (b + 2*c*theta_f + 3*d*theta_f^2); 
            mju/(r2_opt^4*theta2_dot_opt^2) - (1/r2_opt + 2*c + 6*d*theta_f)];

    efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
    
    a_n = a;
    b_n = b;
    c_n = c;
    e_n = double(subs(efg(1), d, d_optimized_tof));
    f_n = double(subs(efg(2), d, d_optimized_tof));
    g_n = double(subs(efg(3), d, d_optimized_tof));
    d_n = double(subs(d, d, d_optimized_tof));
    
    theta_n = double(subs(theta, theta, linspace(theta_0, theta_f, plotAccuracy)));
    r_es = 1 ./ (a_n + b_n*theta_n + c_n*theta_n.^2 + d_n*theta_n.^3 + e_n*theta_n.^4 + f_n*theta_n.^5 + g_n*theta_n.^6);
    x = cos(theta_n+theta1) .* r_es;
    y = sin(theta_n+theta1) .* r_es;
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    r_n = 1 / (a_n + b_n*theta + c_n*theta^2 + d_n*theta^3 + e_n*theta^4 + f_n*theta^5 + g_n*theta^6);
    timeFunction_opt = sqrt((r_n^4/mju) * (1/r_n + 2*c_n + 6*d_n*theta + 12*e_n*theta^2 + 20*f_n*theta^3 + 30*g_n*theta^4));

    timeFunction_opt_nn = @(angle) double(subs(timeFunction_opt, theta, angle));

    dv_optimized_tof = trapz(theta_n, timeFunction_opt_nn(theta_n));

    title(sprintf("deltaV optimized trajectory\n%.0f m/s in %s ", deltaResult, secToTime(dv_optimized_tof)));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end

%% Plotting thrust curves
figure;
hold on;
if optimizeDV == 1
    % Plot results of dV optimized trajectory
    plot(thrustCurve_dV(1, :), thrustCurve_dV(2, :));
end
if optimizeTOF == 1
    % Plot results of dV optimized trajectory
    plot(thrustCurve_tof(1, :), thrustCurve_dV(2, :));
end
  
title(sprintf("Comparision of thrust curves for TOF solution and dV optimization"));
legend("deltaV Optimized trajectory", "TOF solution trajectory");  
xlabel("theta");
ylabel("thurst [N]");


%% Second to time conversion
function [timestring] = secToTime(seconds)
    years = floor(seconds / (365.25 * 86400));
    seconds = seconds - years * (365.25 * 86400);
    days = floor(seconds / 86400);
    seconds = seconds - days * 86400;
    hours = floor(seconds / 3600);
    seconds = seconds - hours * 3600;
    mins = floor(seconds / 60);
    seconds = seconds - mins * 60;
    secs = floor(seconds);

    timestring = sprintf("%.0f y %.0f d %.0f h %.0f m %.0f s", years, days, hours, mins, secs);
end
