%% Shape-Based Approach to Low-Thrust Rendezvous Trajectory Design
clear all;
close all;

%Things to change to test stuff out:
% tf, desired time of flight
% N, number of rotations around central body
% theta_tilde, angle between initial and end position
% Any of the orbital parameters
% limits for d
% And stuff

%% Declaring globals

global theta_f theta_0 intApprox a_initial a_final;
global currentTime transferDateGuess mju N tof_current;
global deltaResult theta2_opt r2_opt theta1_opt r1_opt;
global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt;
global theta2 omega1 omega2 e1 e2 theta1 r1 r2;
global thrustFunction thetaDotFunction radiusFunction;
global d_solution timeResult dateOptimal ;
global Tp1 Tp2 TOF_estimation d_opt timeFunction_nn;
global thrustFunction_opt thetaDotFunction_opt%; thetaDotSquareFunction_opt;
global timeFunction_opt radiusFunction_opt plotAccuracy;
global d_minimum d_maximum rMin rMax theta_vec initial_DeltaV;
global nu2_i r2_i

global plotDV_3D

%% Estimation values and set parameters

%--First Orbit Parameters--
%Semimajor axis
%a_initial= 6800*1000;
a_initial = 200*10^9;
%Time of last perigee pass
Tp1 = 0;
%Eccentricity
e1 = 0.5;
%Argument of perigee
omega1 = 0;

%--Second Orbit Parameters--
%Semimajor axis
%a_final = 42164*1000;
a_final = 500*10^9;
%Time of last perigee pass for object 2
Tp2 = 0;
%Eccentricity
e2 = 0.7;
%Argument of perigee
omega2 = pi;

%Gravitational parameter
%mju = 3.986004418*10^14; %Earth
mju = 1.32712440018*10^20; %Sun

%Minimum allowed radius
%rMin = (6371 + 250) * 1e3; %Earth
rMin = (696340 + 100000) * 1e3; %Sun

%Maximum allowed radius
rMax = 5*max(a_initial, a_final);

%Number of rotations around central body
N = 1;

%Solve TOF, optimize deltaV and/or optimize transfer date
optimizeTOF = 1;
optimizeDV = 1;
optimizeDATE = 1;

%Accuracies of approximation
intApprox = 50;
plotAccuracy = 1000;

%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

% CHANGED - Target will be at random point along 2nd orbit
initialTime = 0;
currentTime = initialTime;   
transferDateGuess = 0;

%Period of the first orbit
P1 = 2*pi/sqrt(mju/a_initial^3);

%How far into the future optimal date is searched.
dateSearchSpan = 3 * P1;

%How many initial points the global search starts with
gsPointCount = 20;

%Symbolic variables
syms theta %Angle along the transfer orbit
syms d %The coefficient used to tune the orbit

%% Calculating the orbital parameters

updateParameters(1);

%Propagating the two orbits
p1 = a_initial * (1-e1^2);
p2 = a_final * (1-e2^2);
nu = linspace(0, 2*pi, plotAccuracy);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

%% TOF optimization + result plotting
if optimizeTOF == 1    
    %Initialize best values
    timeResult = Inf;

    opt = optimset('TolFun', 1e-18, 'TolX', 1e-18);

    %TOF solution
    d_solution = fzero(@transferTimeOptimization, [d_minimum, d_maximum], opt);

    %Plot results of TOF solved trajectory
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    % TOF solving does not update parameters --> No need to use the saved
    % optimal functions and values
    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_solution);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    
    thrustFunction_n = subs(thrustFunction, d, d_solution);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));

    radiusFunction_n = subs(radiusFunction, d, d_solution);
    radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));
    
    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);

    thrustCurve_tof = [theta_vec_plot; thrustFunction_nn(theta_vec_plot)];    

    deltaV_tof = trapz(theta_vec_plot, abs(jerkFunction_nn(theta_vec_plot)));
    initial_DeltaV = deltaV_tof;

    x = cos(theta_vec_plot+theta1) .* radiusFunction_nn(theta_vec_plot);
    y = sin(theta_vec_plot+theta1) .* radiusFunction_nn(theta_vec_plot);
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    title(sprintf("TOF optimized trajectory for fixed starting point\nTarget TOF: %s\nAchieved TOF: %s\nRequired deltaV: %.0f m/s", secToTime(TOF_estimation), secToTime(TOF_estimation+timeResult), deltaV_tof));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
    
end
%% Delta V optimization for fixed starting point + result plotting
if optimizeDV == 1
    plotDV_3D = 0;
    %Initialize best values
    deltaResult = Inf;
    timeResult = Inf;
    
    %Reset current Time
    currentTime = initialTime;
    tof_current = TOF_estimation;

    %Trying out process plotting:
    figure;
    hold on;

    opt = optimset('TolFun',1e2, 'TolX', 1e3);
    tf_fuelOptimal = fminsearch(@optimalDVSolver, TOF_estimation, opt);

    jerkFunction_n = subs(thrustFunction_opt/thetaDotFunction_opt, d, d_opt);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    thrustFunction_n = subs(thrustFunction_opt, d, d_opt);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));
    timeFunction_n = subs(timeFunction_opt, d, d_opt);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    radiusFunction_n = subs(radiusFunction_opt, d, d_opt);
    radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));
    
    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);
    
    thrustCurve_dV = [theta_vec_plot; thrustFunction_nn(theta_vec_plot)]; 
    
    deltaV_opt = trapz(theta_vec_plot, abs(jerkFunction_nn(theta_vec_plot)));
    initial_DeltaV = deltaV_opt;

    time_t = trapz(theta_vec_plot, timeFunction_nn(theta_vec_plot));
    
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega1 + nu1_i_opt) * r1_i_opt, sin(omega1 + nu1_i_opt) * r1_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','c')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* radiusFunction_nn(theta_vec_plot);
    y = sin(theta_vec_plot+theta1_opt) .* radiusFunction_nn(theta_vec_plot);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized trajectory\nAchieved TOF: %s\n%.0f : %.0f m/s", secToTime(time_t), deltaV_opt, deltaResult));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end
%% Optimize the transfer date for deltaV + result plotting
if optimizeDATE == 1
    plotDV_3D = 1;

    %Initialize best values
    deltaResult = Inf;
    timeResult = Inf;
    
    %Reset current Time
    tof_current = TOF_estimation;

    figure;
    hold on;


    for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
        currentTime = time;
       
        %Try plotting optimization

%         gs = MultiStart;
%         numberOfStartPoints = 20;
%         opt = optimset('TolFun',1e2, 'TolX', 1e5);
%         problem = createOptimProblem('fmincon','x0',[initialTime, TOF_estimation],...
%         'objective',@transferDateOptimization,'lb',[initialTime, 0.1*TOF_estimation], ...
%         'ub',[initialTime + 3*P1, 10*TOF_estimation], 'options', opt);
%         x = run(gs,problem, numberOfStartPoints);
    
        opt = optimset('TolFun',1e2, 'TolX', 1e5);
        tf_fuelOptimal = fminsearch(@optimalDVSolver, TOF_estimation, opt);
    end

    jerkFunction_n = subs(thrustFunction_opt/thetaDotFunction_opt, d, d_opt);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    thrustFunction_n = subs(thrustFunction_opt, d, d_opt);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));
    timeFunction_n = subs(timeFunction_opt, d, d_opt);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    radiusFunction_n = subs(radiusFunction_opt, d, d_opt);
    radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_date = [theta_vec_plot; thrustFunction_nn(theta_vec_plot)]; 

    deltaV_opt = trapz(theta_vec_plot, abs(jerkFunction_nn(theta_vec_plot)));
    time_t = trapz(theta_vec_plot, timeFunction_nn(theta_vec_plot));

    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega1 + nu1_i_opt) * r1_i_opt, sin(omega1 + nu1_i_opt) * r1_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','c')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* radiusFunction_nn(theta_vec_plot);
    y = sin(theta_vec_plot+theta1_opt) .* radiusFunction_nn(theta_vec_plot);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized transfer date\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s",secToTime(dateOptimal), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end

%% Plotting thrust curves
figure;
hold on;
if optimizeTOF == 1
    % Plot results of tof optimized trajectory
    plot(thrustCurve_tof(1, :), thrustCurve_tof(2, :));
    d_v_2 = trapz(thrustCurve_tof(1, :), abs(thrustCurve_tof(2, :)));
end
if optimizeDV == 1
    % Plot results of dV optimized trajectory
    plot(thrustCurve_dV(1, :), thrustCurve_dV(2, :));
    d_v_1 = trapz(thrustCurve_dV(1, :), abs(thrustCurve_dV(2, :)));
end
if optimizeDATE == 1
    % Plot results of date optimized trajectory
    plot(thrustCurve_date(1, :), thrustCurve_tof(2, :));
    d_v_2 = trapz(thrustCurve_tof(1, :), abs(thrustCurve_tof(2, :)));
end

title(sprintf("Comparision of thrust curves for all solutions"));
legend("TOF solution trajectory", "deltaV Optimized TOF", "deltaV Optimized TOF and Date");  
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
