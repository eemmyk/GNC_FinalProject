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

%% Declaring globals

global theta_0 intApprox a_initial a_final;
global currentTime previousTime mju N tof_current;
global deltaResult;
global omega1 omega2 e1 e2;
global d_solution dateOptimal tof_optimal;
global Tp1 Tp2 TOF_estimation d_opt;
global plotAccuracy;
global d_minimum d_maximum rMin rMax initial_DeltaV;
global n1 P1 n2 P2 p1 p2

global safeTransferAngle TOF_corrMult dAdjustment;

global solveDate plotTransferWindow

global opt_nu_fzero opt_tf_angle opt_d_lim_fzero opt_tof_fzero;

global theta_vec

global paramVector paramVector_opt;

global tfWindowPixelsX tfWindowPixelsY;

% global contourMap contIndX contIndY contIndLim


%% paramVector(_opt) structure:
% 1:  mju
% 2:  gamma1(_opt)
% 3:  gamma2(_opt)
% 4:  theta_f(_opt)
% 5:  theta1_dot(_opt)
% 6:  theta2_dot(_opt)
% 7:  r1(_opt)
% 8:  r2(_opt)
% 9:  theta1(_opt)
% 10: theta2(_opt)
% 11: nu2_i(_opt)
% 12: r2_i(_opt)

%% Central body information

%Gravitational parameter
%mju = 3.986004418*10^14; %Earth
mju = 1.32712440018*10^20; %Sun

%How large is the central body
%bodyRadius = 6371 * 1e3; %[km] Earth
bodyRadius = 696340 * 1e3; %[km] Sun

%Minimum allowed radius
%rMin = bodyRadius + 250 * 1e3; %[km] Earth
rMin = bodyRadius + 100000 * 1e3; %[km] Sun

%% Estimation values and set parameters

%--First Orbit Parameters--
%Semimajor axis
%a_initial= 6800*1000;
a_initial = 225*10^9;
%Period of the orbit
P1 = 2*pi/sqrt(mju/a_initial^3);
%Time of last perigee pass
Tp1 = 0;
%Eccentricity
e1 = 0;
%Argument of perigee
omega1 = 0;

%--Second Orbit Parameters--
%Semimajor axis
%a_final = 384748*1000;
a_final = 150*10^9;
%Period of the orbit
P2 = 2*pi/sqrt(mju/a_final^3);
%Time of last perigee pass
Tp2 = 0;
%Eccentricity
e2 = 0;
%Argument of perigee
omega2 = 0;

%Number of additional rotations around central body
N = 0;
  
%Limits for TOF in relation to first guess
TofLimLow = 0.3;
TofLimHigh = 5;

%--Adjustment Parameters--
%How much the calculated TOF is increased
TOF_corrMult = 1;
%How much d-coefficients are changed to find new solutions
dAdjustment = 1.2;

%Minimum angle travelled around the central body
safeTransferAngle = 0;%pi/2;

%Maximum allowed radius
rMax = 10*max(a_initial, a_final);

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

%Testing out lcm to always find interesting values
years1 = ceil(P1/(365*86400));
years2 = ceil(P2/(365*86400));
dateSearchSpan = lcm(years1, years2) * 365 * 86400;

%How many initial points the global search starts with
%Linear for option 1
%Linear for option 2
%Squared for option 3
gsPointCount = 64;

%Which approach to global search is taken
%Option 1: Global search
%Option 2: TOF search with date vector
%Option 3: TOF and date vectors
transferWindowSearchOption = 3;

%Set up the initial deltaV value for plotting
initial_DeltaV = 1e24; %A big number

% contourMap = zeros(gsPointCount);
% contIndX = 1;
% contIndY = 1;
% contIndLim = gsPointCount;

%% Define some optimization options here for speeeeeeeeed!
opt_tof_fzero = optimset('TolFun', 1e2, 'Display', 'off');
opt_tof_fzero_acc = optimset('TolFun', 1e-1);
opt_nu_fzero = optimset('TolFun', 1e-1, 'TolX', 1e-1, 'Display', 'off');
opt_tf_angle = optimset('TolFun', 1e-3, 'Display', 'off');
opt_d_lim_fzero = optimset('TolFun', 1e2, 'TolX', 1e-15, 'Display', 'off');
opt_dv_fminsearch = optimset('TolFun',1e2, 'TolX', min(P1, P2)/1000);
opt_dv_global = optimset('TolFun',1e2, 'TolX', 1e4);

%% Calculating the orbital parameters
n1 = sqrt(mju/a_initial^3);
p1 = a_initial * (1-e1^2);

n2 = sqrt(mju/a_final^3);
p2 = a_final * (1-e2^2);

updateParameters(1);

%How large are the pixels in the transfer window plot
tfWindowPixelsX = ceil(dateSearchSpan / (gsPointCount-1));
tfWindowPixelsY = ceil((TofLimHigh - TofLimLow) * TOF_estimation / (gsPointCount-1));

theta_f = paramVector(4);
r1 = paramVector(7);
r2 = paramVector(8);
theta1 = paramVector(9);
theta2 = paramVector(10);
nu2_i = paramVector(11);
r2_i = paramVector(12);

%Propagating the two orbits
nu = linspace(0, 2*pi, plotAccuracy);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

%% TOF optimization + result plotting
if optimizeTOF == 1    
    try
        %TOF solution
        tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, tof_current, theta_vec);
        d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], opt_tof_fzero_acc);
    catch
        fprintf("Inital Time of Flight guess not achievable\n")
        d_solution = (d_minimum + d_maximum)/2;
    end
    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);

    %Zero at the end because no optimize values exist yet
    thrustCurve_tof = [theta_vec_plot; fThrustFunction(d_solution, theta_vec_plot, paramVector)];    

    deltaV_tof = trapz(theta_vec_plot, abs(fJerkFunction(d_solution, theta_vec_plot, paramVector)));
    initial_DeltaV = deltaV_tof;
    
    time_t = trapz(theta_vec_plot, fTimeFunction(d_solution, theta_vec_plot, paramVector));

    %Plot results of TOF solved trajectory
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
        
    rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, paramVector);
    y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, paramVector);
    
    plot(x, y, "Color", [0.2 0.7 0.2]);
    
    title(sprintf("TOF solution trajectory\nTarget TOF: %s\nAchieved TOF: %s\nRequired deltaV: %.0f m/s", secToTime(TOF_estimation), secToTime(time_t), deltaV_tof));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal

end
%% Delta V optimization for fixed starting point + result plotting
if optimizeDV == 1
    solveDate = 0;
    plotTransferWindow = 0;

    %Initialize best values
    deltaResult = Inf;

    tf_fuelOptimal = fminsearch(@optimalDVSolver, TOF_estimation, opt_dv_fminsearch);

    %Update local variables to the optimal solution
    theta_f_opt = paramVector_opt(4);
    r1_opt = paramVector_opt(7);
    r2_opt = paramVector_opt(8);
    theta1_opt = paramVector_opt(9);
    theta2_opt = paramVector_opt(10);
    nu2_i_opt = paramVector_opt(11);
    r2_i_opt = paramVector_opt(12);

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_dV = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, paramVector_opt)]; 
    
    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, paramVector_opt)));
    initial_DeltaV = deltaV_opt;

    time_t = trapz(theta_vec_plot, fTimeFunction(d_opt, theta_vec_plot, paramVector_opt));
    
    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
    
    bodyScale = max(a_initial, a_final) * 0.075;
    
    rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, paramVector_opt);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, paramVector_opt);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("deltaV optimized trajectory\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s", secToTime(initialTime), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end
%% Optimize the transfer date for deltaV + result plotting
if optimizeDATE == 1
    plotTransferWindow = 1;
%     tw_graph_ind = 1;
%     tw_graph = zeros(gsPointCount*gsPointCount,6);

    %Initialize best values
    deltaResult = Inf;
    
    %Reset current Time
    tof_current = TOF_estimation;

    if plotTransferWindow == 1
        figure;
        hold on;
    end

    %-- Solve transfer window using Global Search --
    if transferWindowSearchOption == 1
        gs = MultiStart;
        numberOfStartPoints = gsPointCount;
        problem = createOptimProblem('fmincon','x0',[initialTime, TOF_estimation],...
        'objective',@optimalDVSolver,'lb',[initialTime, TofLimLow*TOF_estimation], ...
        'ub',[initialTime + dateSearchSpan, TofLimHigh*TOF_estimation], 'options', opt_dv_global);
        run(gs, problem, numberOfStartPoints);
    end    
    
    %-- Solve transfer window using set dates and optimize TOF --
    if transferWindowSearchOption == 2
        for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong> \n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            %Try plotting optimization
            solveDate = 0;
            currentTime = time;
            fminsearch(@optimalDVSolver, TOF_estimation, opt_dv_fminsearch);
            %tof_current = tof_optimal;
        end
    end
    %-- Solve transfer window using a set grid of dates and TOFs --
    if transferWindowSearchOption == 3
        solveDate = 1;
        for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong> \n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof = linspace(TofLimLow*TOF_estimation, TofLimHigh*TOF_estimation, gsPointCount)
    
                currentTime = time;
                tof_current = tof;
        
                optimalDVSolver([currentTime, tof_current]);
            end
        end
    end
    pause(0)
    
    %All searches are complimented by a final search for the local minimum 
    %close to the best found solution
    solveDate = 1;
    plotTransferWindow = 0;
    tf_fuelOptimal = fminsearch(@optimalDVSolver, [dateOptimal, tof_optimal], opt_dv_fminsearch);
    fprintf("Optimization Completed! ---- Finalized %cV Found: <strong>%.0f m/s</strong> \n", 916, deltaResult);

    %contourf(contourMap, 100);
    if plotTransferWindow == 1
        %plot3(tw_graph(:,1), tw_graph(:,2), tw_graph(:,3),'o','Color','k','MarkerSize',10)
        title(sprintf("DeltaV Map For Different Launch Dates and Time of Flights"));
        xlabel("Time past initial Time [s]");
        ylabel("Time of Flight [s]")
    end

    %Update local variables to the optimal solution
    theta_f_opt = paramVector_opt(4);
    r1_opt = paramVector_opt(7);
    r2_opt = paramVector_opt(8);
    theta1_opt = paramVector_opt(9);
    theta2_opt = paramVector_opt(10);
    nu2_i_opt = paramVector_opt(11);
    r2_i_opt = paramVector_opt(12);

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_date = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, paramVector_opt)]; 

    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, paramVector_opt)));

    time_t = trapz(theta_vec_plot, fTimeFunction(d_opt, theta_vec_plot, paramVector_opt));

    figure;
    hold on;
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
        
    rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, paramVector_opt);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_plot, paramVector_opt);
   
    plot(x, y, "Color", [0.2 0.7 0.2]);

    title(sprintf("Full transfer window solution\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s",secToTime(dateOptimal), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
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
