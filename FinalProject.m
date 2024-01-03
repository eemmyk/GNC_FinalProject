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
global currentTime previousTime mju N tof_current;
global deltaResult theta2_opt r2_opt theta1_opt r1_opt;
global nu2_i_opt r2_i_opt theta_f_opt;
global theta2 omega1 omega2 e1 e2 theta1 r1 r2;
global d_solution dateOptimal ;
global Tp1 Tp2 TOF_estimation d_opt;
global plotAccuracy;
global d_minimum d_maximum rMin rMax initial_DeltaV;
global nu2_i r2_i;
global n1 P1 n2 P2 p1 p2

global safeTransferAngle TOF_corrMult;

global solveDate plot3D

% global contourMap contIndX contIndY contIndLim

%% Central body information

%Gravitational parameter
%mju = 3.986004418*10^14; %Earth
mju = 1.32712440018*10^20; %Sun

%Minimum allowed radius
%rMin = (6371 + 250) * 1e3; %Earth
rMin = (696340 + 100000) * 1e3; %Sun

%% Estimation values and set parameters

%--First Orbit Parameters--
%Semimajor axis
%a_initial= 6800*1000;
a_initial = 150*10^9;
%Period of the orbit
P1 = 2*pi/sqrt(mju/a_initial^3);
%Time of last perigee pass
Tp1 = 0;
%Eccentricity
e1 = 0.5;
%Argument of perigee
omega1 = -pi/3;

%--Second Orbit Parameters--
%Semimajor axis
%a_final = 42164*1000;
a_final = 225*10^9;
%Period of the orbit
P2 = 2*pi/sqrt(mju/a_final^3);
%Time of last perigee pass for object 2
Tp2 = 0;
%Eccentricity
e2 = 0.6;
%Argument of perigee
omega2 = -pi/2;

%Limits for TOF in relation to first guess
TofLimLow = 0.2;
TofLimHigh = 5;
%How much the calculated TOF is increased
TOF_corrMult = 1;

%Minimum angle travelled around the central body
safeTransferAngle = pi/2;

%Maximum allowed radius
rMax = 10*max(a_initial, a_final);

%Number of rotations around central body
N = 0;

%Spacecraft mass [16U CubeSat]
m = 32; %kg

%Solve TOF, optimize deltaV and/or optimize transfer date
optimizeTOF = 1;
optimizeDV = 1;
optimizeDATE = 1;

%Accuracies of approximation
intApprox = 100;
plotAccuracy = 1000;

%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

% CHANGED - Target will be at random point along 2nd orbit
initialTime = 0;
currentTime = initialTime;
%Impossible time to trigger parameter generation on first run
%if previousTime ~= currentTime
previousTime = -1;

%How far into the future optimal date is searched.
dateSearchSpan = 2*P1;%4*max(P1, P2);

%How many initial points the global search starts with
gsPointCount = 48;

% contourMap = zeros(gsPointCount);
% contIndX = 1;
% contIndY = 1;
% contIndLim = gsPointCount;


%% Calculating the orbital parameters
n1 = sqrt(mju/a_initial^3);
p1 = a_initial * (1-e1^2);

n2 = sqrt(mju/a_final^3);
p2 = a_final * (1-e2^2);

updateParameters(1);

%Propagating the two orbits
nu = linspace(0, 2*pi, plotAccuracy);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

%% TOF optimization + result plotting
if optimizeTOF == 1    
    opt = optimset('TolFun', 1e-1);
    
    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);

    time_t_max = trapz(theta_vec_plot, fTimeFunction(d_minimum, theta_vec_plot, 0));
    time_t_min = trapz(theta_vec_plot, fTimeFunction(d_maximum, theta_vec_plot, 0));
    
%     figure;
%     hold on;
%     plot(real(fTimeFunction(d_minimum, theta_vec_plot, 0)));
%     plot(real(fTimeFunction(d_maximum, theta_vec_plot, 0)));

    %TOF solution
    d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt);

    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);

    %Zero at the end because no optimize values exist yet
    thrustCurve_tof = [theta_vec_plot; fThrustFunction(d_solution, theta_vec_plot, 0)];    

    deltaV_tof = trapz(theta_vec_plot, abs(fJerkFunction(d_solution, theta_vec_plot, 0)));
    %title("Gamma-TOF")

    initial_DeltaV = deltaV_tof;
    
    time_t = trapz(theta_vec_plot, fTimeFunction(d_solution, theta_vec_plot, 0));

    x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
    y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
    
    %figure;
%     for time = linspace(initialTime, initialTime + dateSearchSpan, 1000)
%         currentTime = time;
%         updateParameters(0)
%         tof_current = TOF_estimation;
%         hold off
%         plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
%         hold on
%         plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
%         
%         plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
%         plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
%         
%         bodyScale = max(a_initial, a_final) * 0.075;
% 
%         rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
%         
%         d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt);
% 
%         theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);
%  
%         x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
%         y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
%         
%         plot(x, y, "Color", [0.2 0.7 0.2]);
% 
%         pause(0)
%     end

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
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    title(sprintf("TOF solution trajectory\nTarget TOF: %s\nAchieved TOF: %s\nRequired deltaV: %.0f m/s", secToTime(TOF_estimation), secToTime(time_t), deltaV_tof));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal

%     figure;
%     diffRad = diff(fRadiusFunction(d_solution, theta_vec_plot, 0));
%     plot(diffRad);
%     title("Diff-TOF");
% 
%     figure;
%     diffDiffRad = diff(diff(fRadiusFunction(d_solution, theta_vec_plot, 0)));
%     diffDiffRad(end+1) = 0;
%     plot(diffDiffRad);
%     title("DiffDiff-TOF");
% 
%     figure;
%     plot(diffRad .* diffDiffRad);
%     title("DiffMult-TOF")

end
%% Delta V optimization for fixed starting point + result plotting
if optimizeDV == 1
    solveDate = 0;
    plot3D = 0;
    %Initialize best values
    deltaResult = Inf;

%     %Reset current Time
%     currentTime = initialTime;
%     tof_current = TOF_estimation;

    %Trying out process plotting:
    figure;
    hold on;

    opt = optimset('TolFun',1e2, 'TolX', min(P1, P2)/1000);
    tf_fuelOptimal = fminsearch(@optimalDVSolver, TOF_estimation, opt);

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_dV = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, 1)]; 
    
    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, 1)));
    %title("Gamma-DV")

    %initial_DeltaV = deltaV_opt;

    time_t = trapz(theta_vec_plot, fTimeFunction(d_opt, theta_vec_plot, 1));
    
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    %plot(cos(omega1 + nu1_i_opt) * r1_i_opt, sin(omega1 + nu1_i_opt) * r1_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','c')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, 1);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, 1);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized trajectory\nLaunch date: %s\nAchieved TOF: %s\n%.0f m/s", secToTime(initialTime), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end
%% Optimize the transfer date for deltaV + result plotting
if optimizeDATE == 1
    solveDate = 1;
    plot3D = 1;


    %Initialize best values
    deltaResult = Inf;
    
    %Reset current Time
    tof_current = TOF_estimation;

    figure;
    hold on;

%     gs = MultiStart;
%     numberOfStartPoints = gsPointCount^2;
%     opt = optimset('TolFun',1e2, 'TolX', 1e4);
%     problem = createOptimProblem('fmincon','x0',[initialTime, TOF_estimation],...
%     'objective',@optimalDVSolver,'lb',[initialTime, TofLimLow*TOF_estimation], ...
%     'ub',[initialTime + dateSearchSpan, TofLimHigh*TOF_estimation], 'options', opt);
%     x = run(gs,problem, numberOfStartPoints);
%     

%     for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
%         %Try plotting optimization
%         solveDate = 0;
%         currentTime = time;
% 
%         opt = optimset('TolFun',1e2, 'TolX', 1e4);
%         tf_fuelOptimal = fminsearch(@optimalDVSolver, tof_current, opt);
% 
%         %tof_current = tof_optimal;
%     end

    

    for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
        fprintf("Optimizing: %.2f %% Completed ---- Best %cV Solution Found: %.0f m/s\n", 100 * time/(initialTime + dateSearchSpan), 916, deltaResult);
        for tof = linspace(TofLimLow*TOF_estimation, TofLimHigh*TOF_estimation, gsPointCount)
        %Try plotting optimization
        currentTime = time;
        tof_current = tof;

        opt = optimset('TolFun',1e2, 'TolX', min(P1, P2)/1000);
        deltaV = optimalDVSolver([currentTime, tof_current]);
        end
    end

    %contourf(contourMap, 100);

    title(sprintf("DeltaV Map For Different Launch Dates and Time of Flights"));
    xlabel("Time past initial Time [s]");
    ylabel("Time of Flight [s]")
    
    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_date = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, 1)]; 

    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, 1)));
    %title("Gamma-DATE")


    time_t = trapz(theta_vec_plot, fTimeFunction(d_opt, theta_vec_plot, 1));

    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    %plot(cos(omega1 + nu1_i_opt) * r1_i_opt, sin(omega1 + nu1_i_opt) * r1_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','c')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, 1);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, 1);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized transfer date\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s",secToTime(dateOptimal), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal


%     figure;
%     diffRad = diff(fRadiusFunction(d_opt, theta_vec_plot, 1));
%     plot(diffRad);
%     title("Diff-DATE");
% 
%     figure;
%     diffDiffRad = diff(diff(fRadiusFunction(d_opt, theta_vec_plot, 1)));
%     diffDiffRad(end+1) = 0;
%     plot(diffDiffRad);
%     title("DiffDiff-DATE");
% 
%     figure;
%     plot(diffRad .* diffDiffRad);
%     title("DiffMult-DATE")
end

%% Plotting thrust curves
figure;
hold on;
if optimizeTOF == 1
    % Plot results of tof optimized trajectory
    plot(thrustCurve_tof(1, :), m*thrustCurve_tof(2, :));
end
if optimizeDV == 1
    % Plot results of dV optimized trajectory
    plot(thrustCurve_dV(1, :), m*thrustCurve_dV(2, :));
end
if optimizeDATE == 1
    % Plot results of date optimized trajectory
    plot(thrustCurve_date(1, :), m*thrustCurve_date(2, :));
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
