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
%global gamma2_opt theta2_dot_opt gamma1_opt theta1_dot_opt;
global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt;
global theta2 omega1 omega2 e1 e2 theta1 r1 r2;
global timeFunction thrustFunction thetaDotFunction thetaDotSquareFunction radiusFunction;
global d_solution timeResult costResult dateOptimal tof_optimal;
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
a_initial = 100*10^9;
%Time of last perigee pass
Tp1 = 0;
%Eccentricity
e1 = 0.2;
%Argument of perigee
omega1 = pi/5;

%--Second Orbit Parameters--
%Semimajor axis
%a_final = 42164*1000;
a_final = 500*10^9;
%Time of last perigee pass for object 2
Tp2 = 0;
%Eccentricity
e2 = 0.1;
%Argument of perigee
omega2 = 3/4*pi;

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

% % 1 + 2N hohmann transfer orbits in time
% TOF_estimation = (1+2*N)*pi*sqrt((a_initial+a_final)^3/(8*mju));
% 
% tof_current = TOF_estimation;

%Used for dV optimization

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
gsPointCount = 10;

%Symbolic variables
syms theta %Angle along the transfer orbit
syms d %The coefficient used to tune the orbit

%% Calculating the orbital parameters

updateParameters(1);

% %initial orbit parameters
% a1 = a_initial;
% 
% T1_i = mod(currentTime - Tp1, P1);
% 
% syms nu_time
% 
% E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e1)/(1-e1)));
% nuSolver = n1*T1_i==pi+E-e1*sin(E);
% nuSolutions_i = vpasolve(nuSolver, nu_time);
% nuSolutions_i = double(nuSolutions_i);
% 
% nu1_i = mod(nuSolutions_i + 2*pi, 2*pi);
% 
% %In reference coords
% theta1 = nu1_i + omega1;
% gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));
% p1 = a1 * (1-e1^2);
% r1 = p1 / (1+e1*cos(nu1_i));
% r1_i = p1 / (1+e1*cos(nu1_i));
% theta1_dot = sqrt(mju/a1^3) * a1^2/r1^2 * sqrt(1-e1^2);
%     
% %Target orbit parameters
% %Some are known, while others are calculated from desired tof and initial
% %conditions
% a2 = a_final;
% 
% 
% %Trying this one instead
% n2 = sqrt(mju/a2^3);
% P2 = 2*pi/n2;
% 
% %Tp2 = rand() * P2;
% 
% T = mod(currentTime - Tp2, P2);
% %T = currentTime - Tp2;
% 
% % nT = E - e*sin(E)
% % E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));
% 
% syms nu_time
% 
% E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
% nuSolver = n2*T==pi+E-e2*sin(E);
% nuSolutions_i = vpasolve(nuSolver, nu_time);
% nuSolutions_i = double(nuSolutions_i);
% 
% nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);
% 
% %Now we can continue by calculating the true anomaly when the spacecraft
% %reaches orbit 2. At currentTime + tf
% 
% finalTime = currentTime + tof_current;
% T = mod(finalTime - Tp2, P2);
% %T = finalTime - Tp2;
% 
% % nT = E - e*sin(E)
% % E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));
% 
% E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
% nuSolver = n2*T==pi+E-e2*sin(E);
% nuSolutions_f = vpasolve(nuSolver, nu_time);
% nuSolutions_f = double(nuSolutions_f);
% 
% nu2 = mod(nuSolutions_f + 2*pi, 2*pi);
% %In reference coords
% theta2 = nu2 + omega2;
% theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
% 
% gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
% p2 = a2 * (1-e2^2);
% r2 = p2 / (1+e2*cos(nu2));
% r2_i = p2 / (1+e2*cos(nu2_i));
% theta2_dot = sqrt(mju/a2^3)*a2^2/r2^2 * sqrt(1-e2^2);
% 
% %The total transfer angle is represented by:
% theta_f = 2*pi * N + theta_tilde;

%Propagating the two orbits
p1 = a_initial * (1-e1^2);
p2 = a_final * (1-e2^2);
nu = linspace(0, 2*pi, plotAccuracy);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

theta_vec = linspace(theta_0, theta_f, intApprox);

%% Solving the coefficients
% syms d theta;
% 
% a = 1/r1;
% b = -tan(gamma1) / r1;
% c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);
% 
% efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
%             -48*theta_f     18*theta_f^2 -2*theta_f^3; 
%              20            -8*theta_f     theta_f^2];
% 
% efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
%             -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
%             mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
% 
% efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
% 
% e = efg(1);
% f = efg(2);
% g = efg(3);
% 
% r = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);
% gamma = atan(-r * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4 + 6*g*theta^5));
% 
% 
% 
% thetaDotFunction = sqrt((mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
% thrustFunction = -mju / (2 * r^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta^2 + 120*g*theta^3 - tan(gamma)/r) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)^2;
% 
% timeFunction = sqrt((r^4/mju) * (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
% timeFunction_nn = @(d_coeff, angle) double(subs(subs(timeFunction, d, d_coeff), theta, angle));
% 
% thetaDotSquareFunction = (mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4);
% 
% radiusFunction = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);

%% Check if TOF is solvable:
% safeTransferAngle = pi;
% 
% if theta_f < safeTransferAngle
%     N = N+1;
%     updateParameters(1);
% 
%     %Used for dV optimization
% end
% %Check what the minimum d coefficient is
% radiusMax_n = subs(1/radiusFunction, theta, theta_f/2);
% radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
% 
% %Move the d-coefficient to the positive direction -->
% %Reduce chance to find orbits going into infinity
% safe_d_multiplier = 1.4;
% d_minimum = fzero(radiusMax_nn, 0);
% if d_minimum < 0
%     d_minimum = d_minimum / safe_d_multiplier;
% else
%     d_minimum = d_minimum * safe_d_multiplier;
% end
% %d_minimum = d_minimum + abs(d_minimum) * safe_d_multiplier;
% 
% %Could additionally be defined as a maximum distance as with the upper
% %limit of the d-coefficient
% 
% radiusMin_n = subs(radiusFunction - rMin, theta, theta_f/2);
% radiusMin_nn = @(d_max_in) double(subs(radiusMin_n, d, d_max_in));
% 
% %Needs to be <1
% d_maximum = fzero(radiusMin_nn, [d_minimum, 1]);

%% TOF optimization + result plotting
if optimizeTOF == 1
    %theta_cross = fminsearch(@thetaSquareFunc, [0, 0]);
    
    %Initialize best values
    timeResult = Inf;

    %currentTime = initialTime;
    %tof_current = TOF_estimation;
    
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
    
    theta_vec = linspace(theta_0, theta_f, plotAccuracy);

    %jerk = [theta_vec; jerkFunction_nn(theta_vec)];    
    thrustCurve_tof = [theta_vec; thrustFunction_nn(theta_vec)];    
    %thetaDotCurve_tof = [theta_vec; thrustFunction_nn(theta_vec) ./ jerkFunction_nn(theta_vec)];

    deltaV_tof = trapz(theta_vec, abs(jerkFunction_nn(theta_vec)));
    initial_DeltaV = deltaV_tof;

    x = cos(theta_vec+theta1) .* radiusFunction_nn(theta_vec);
    y = sin(theta_vec+theta1) .* radiusFunction_nn(theta_vec);
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    title(sprintf("TOF optimized trajectory for fixed starting point\nTarget TOF: %s\nAchieved TOF: %s\nRequired deltaV: %.0f m/s", secToTime(TOF_estimation), secToTime(TOF_estimation+timeResult), deltaV_tof));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
    
end
%end
%% Delta V optimization for fixed starting point + result plotting
if optimizeDV == 1
    plotDV_3D = 0;
    %Initialize best values
    deltaResult = Inf;
    timeResult = Inf;
    %interDeltaResult = Inf;
    
    %Reset current Time
    currentTime = initialTime;
    tof_current = TOF_estimation;
    tf_guess = TOF_estimation;

    %Trying out process plotting:
    figure;
    hold on;

    opt = optimset('TolFun',1e2, 'TolX', 1e5);
    %tf_fuelOptimal = fmincon(@optimalDVSolver, 0, [], [], [], [], 0, [], [], opt);

    tf_fuelOptimal = fminsearch(@optimalDVSolver, tf_guess, opt);

    jerkFunction_n = subs(thrustFunction_opt/thetaDotFunction_opt, d, d_opt);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    thrustFunction_n = subs(thrustFunction_opt, d, d_opt);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));
    timeFunction_n = subs(timeFunction_opt, d, d_opt);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    radiusFunction_n = subs(radiusFunction_opt, d, d_opt);
    radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));
    
    theta_vec = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_dV = [theta_vec; thrustFunction_nn(theta_vec)]; 
    
    deltaV_opt = trapz(theta_vec, abs(jerkFunction_nn(theta_vec)));
    initial_DeltaV = deltaV_opt;

    time_t = trapz(theta_vec, timeFunction_nn(theta_vec));
    
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
    
    x = cos(theta_vec+theta1_opt) .* radiusFunction_nn(theta_vec);
    y = sin(theta_vec+theta1_opt) .* radiusFunction_nn(theta_vec);
   
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
    costResult = Inf;
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
%         %opt = optimset('DiffMinChange', 1e3, 'DiffMaxChange', 1e6);
%         problem = createOptimProblem('fmincon','x0',[initialTime, TOF_estimation],...
%         'objective',@transferDateOptimization,'lb',[initialTime, 0.1*TOF_estimation], ...
%         'ub',[initialTime + 3*P1, 10*TOF_estimation]);%, 'options', opt);
%         x = run(gs,problem, numberOfStartPoints);
    
        opt = optimset('TolFun',1e2, 'TolX', 1e5);
        tf_fuelOptimal = fminsearch(@optimalDVSolver, tf_guess, opt);
    end

    jerkFunction_n = subs(thrustFunction_opt/thetaDotFunction_opt, d, d_opt);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    thrustFunction_n = subs(thrustFunction_opt, d, d_opt);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));
    timeFunction_n = subs(timeFunction_opt, d, d_opt);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
    radiusFunction_n = subs(radiusFunction_opt, d, d_opt);
    radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));

    theta_vec = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_dV = [theta_vec; thrustFunction_nn(theta_vec)]; 

    deltaV_opt = trapz(theta_vec, abs(jerkFunction_nn(theta_vec)));
    time_t = trapz(theta_vec, timeFunction_nn(theta_vec));

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
    
    x = cos(theta_vec+theta1_opt) .* radiusFunction_nn(theta_vec);
    y = sin(theta_vec+theta1_opt) .* radiusFunction_nn(theta_vec);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized transfer date\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s",secToTime(dateOptimal), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end

%% Plotting thrust curves
figure;
hold on;
if optimizeDV == 1
    % Plot results of dV optimized trajectory
    plot(thrustCurve_dV(1, :), thrustCurve_dV(2, :));
    d_v_1 = trapz(thrustCurve_dV(1, :), abs(thrustCurve_dV(2, :)));
end
if optimizeTOF == 1
    % Plot results of tof optimized trajectory
    plot(thrustCurve_tof(1, :), thrustCurve_tof(2, :));
    d_v_2 = trapz(thrustCurve_tof(1, :), abs(thrustCurve_tof(2, :)));
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
