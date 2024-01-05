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
global d_solution dateOptimal tof_optimal;
global Tp1 Tp2 TOF_estimation d_opt;
global plotAccuracy;
global d_minimum d_maximum rMin rMax initial_DeltaV;
global nu2_i r2_i;
global n1 P1 n2 P2 p1 p2

global safeTransferAngle TOF_corrMult dAdjustment;

global solveDate plotTransferWindow %tw_graph_ind tw_graph

global opt_nu_fzero opt_tf_angle opt_d_lim_fzero opt_tof_fzero;

global theta_vec


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
e1 = 0.3;
%Argument of perigee
omega1 = -pi*1.5;

%--Second Orbit Parameters--
%Semimajor axis
%a_final = 42164*1000;
a_final = 350*10^9;
%Period of the orbit
P2 = 2*pi/sqrt(mju/a_final^3);
%Time of last perigee pass for object 2
Tp2 = 0;
%Eccentricity
e2 = 0.6;
%Argument of perigee
omega2 = 0.2*pi;

%Limits for TOF in relation to first guess
TofLimLow = 0.3;
TofLimHigh = 3;

%--Adjustment Parameters--
%How much the calculated TOF is increased
TOF_corrMult = 1.5;
%How much d-coefficients are changed to find new solutions
dAdjustment = 1.1;


%Minimum angle travelled around the central body
safeTransferAngle = pi;

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
intApprox = 50;
plotAccuracy = 1000;

%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

% CHANGED - Target will be at random point along 2nd orbit
initialTime = 45e6;
currentTime = initialTime;
%Impossible time to trigger parameter generation on first run
%if previousTime ~= currentTime
previousTime = -1;

%How far into the future optimal date is searched.
dateSearchSpan = 2*min(P1, P2);

%How many initial points the global search starts with
%Linear for option 1
%Linear for option 2
%Squared (16/9 ratio) for option 3
gsPointCount = 16;

%Which approach to global search is taken
%Option 1: Global search
%Option 2: TOF search with date vector
%Option 3: TOF and date vectors
transferWindowSearchOption = 2;


% contourMap = zeros(gsPointCount);
% contIndX = 1;
% contIndY = 1;
% contIndLim = gsPointCount;

%% Define some optimization options here for speeeeeeeeed!
opt_tof_fzero = optimset('TolFun', 1e2, 'Display', 'off');
opt_tof_fzero_acc = optimset('TolFun', 1e-1);
opt_nu_fzero = optimset('TolFun', 1e-3);
opt_tf_angle = optimset('TolFun', 1e-3, 'Display', 'off');
opt_d_lim_fzero = optimset('TolFun', 1e1, 'TolX', 1e-15, 'Display', 'off');
opt_dv_fminsearch = optimset('TolFun',1e2, 'TolX', min(P1, P2)/1000);
opt_dv_global = optimset('TolFun',1e2, 'TolX', 1e4);

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
    
%     theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);
%     %theta_vec = linspace(theta_0, theta_f, intApprox);
% 
% 
%     time_t_max = trapz(theta_vec_plot, fTimeFunction(d_minimum, theta_vec_plot, 0));
%     time_t_min = trapz(theta_vec_plot, fTimeFunction(d_maximum, theta_vec_plot, 0));
% 
% %     figure;
% %     hold on;
% %     plot(real(fTimeFunction(d_minimum, theta_vec_plot, 0)));
% %     plot(real(fTimeFunction(d_maximum, theta_vec_plot, 0)));
% 
%     figure;
% 
%     solveDate = 0;
%     initial_DeltaV = 10000;
%     deltaResult = Inf;
%     plotTransferWindow = 0;
%     
%     minTof = Inf;
%     maxTof = 0;
% 
%     imaginaryCount = 0;
%     timespanCount = 0;
%     infiniteCount = 0;
%     unNumericCount = 0;
%     unfeasibleCount = 0;
% 
%     for time = linspace(0, (P1+P2), 1000)
%             currentTime = time;
%             updateParameters(1)
%             %tf_fuelOptimal = fminsearch(@optimalDVSolver, tof_current, opt_dv_fminsearch);
% 
% %             deltaV = optimalDVSolver([currentTime, tof]);
%             try 
%                 d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt_tof_fzero_acc);
%                 tfColor = 'green';
%             catch
%                 d_solution = 0;
%                 t_error_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0)) - tof_current
%                 t_error_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0)) - tof_current
%     
%                 imaginarySolutions = ~(isreal(t_error_min) && isreal(t_error_max));
%                 noCrossingSolutions = (t_error_min < 0) == (t_error_max < 0);
%                 infiniteSolutions = ~(isfinite(t_error_min) && isfinite(t_error_max));
%                 unNumericSolutions = isnan(t_error_min) || isnan(t_error_max);
% 
%                 infiniteCount = infiniteCount + infiniteSolutions;
%                 unNumericCount = unNumericCount + unNumericSolutions;
%                 imaginaryCount = imaginaryCount + imaginarySolutions;
%                 %timespanCount = timespanCount + noCrossingSolutions;
% 
%                 if ~imaginarySolutions
%                     timespanCount = timespanCount + noCrossingSolutions;
%                 end
% 
% %                 if ~imaginarySolutions
% %                     tof_current
% %                     t_error_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0))
% %                     t_error_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0))
% %                     tmin = transferTimeSolution(d_maximum)
% %                     tmax = transferTimeSolution(d_minimum)
% %                     
% %                     d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt_tof_fzero_acc);
% %                 end
% 
%                 unfeasibleCount = unfeasibleCount + 1;
%                 tfColor = 'red';
%             end
%             hold off
%             plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
%             hold on
%             plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
%             
%             plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
%             plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
%             
%             bodyScale = max(a_initial, a_final) * 0.075;
%     
%             rectangle('Position',[-0.5*bodyScale, -0.5*bodyScale, bodyScale, bodyScale],'Curvature',[1 1], 'FaceColor',"yellow")
%                 
%             %if deltaResult < 1e23
% 
%                 minTof = min(minTof, tof_current);
%                 maxTof = max(maxTof, tof_current);
%                 
%                 theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);
%                 x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
%                 y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
%                 
%                 plot(x, y, "Color", tfColor);
%             %end
%         pause(0.001)
%     end
% 
%     unfeasibleCount
%     imaginaryCount
%     timespanCount
%     infiniteCount
%     unNumericCount
%     title(sprintf("Minimum possible TOF: %s\nMaximum possible TOF: %s", secToTime(minTof), secToTime(maxTof)));

    %TOF solution
    d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt_tof_fzero_acc);

    theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);

    %Zero at the end because no optimize values exist yet
    thrustCurve_tof = [theta_vec_plot; fThrustFunction(d_solution, theta_vec_plot, 0)];    

    deltaV_tof = trapz(theta_vec_plot, abs(fJerkFunction(d_solution, theta_vec_plot, 0)));

    initial_DeltaV = deltaV_tof;
    
    time_t = trapz(theta_vec_plot, fTimeFunction(d_solution, theta_vec_plot, 0));

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
    
    x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
    y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_plot, 0);
    
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

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_dV = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, 1)]; 
    
    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, 1)));

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

    title(sprintf("deltaV optimized trajectory\nTransfer date: %s\nAchieved TOF: %s\n%.0f m/s", secToTime(initialTime), secToTime(time_t), deltaV_opt));
    legend("Initial orbit", "Target orbit", "body 1 @ t = 0", "body 2 @ t = tf", " body 2 @ t = 0", "Transfer Orbits");
    axis equal
end
%% Optimize the transfer date for deltaV + result plotting
if optimizeDATE == 1
    solveDate = 1;
    plotTransferWindow = 1;
%     tw_graph_ind = 1;
%     tw_graph = zeros(gsPointCount*gsPointCount,6);

    %Initialize best values
    deltaResult = Inf;
    
    %Reset current Time
    tof_current = TOF_estimation;

    figure;
    hold on;

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

            tf_fuelOptimal = fminsearch(@optimalDVSolver, tof_current, opt_dv_fminsearch);
    
            %tof_current = tof_optimal;
        end
    end

    %-- Solve transfer window using a set grid of dates and TOFs --
    if transferWindowSearchOption == 3
        for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong> \n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof = linspace(TofLimLow*TOF_estimation, TofLimHigh*TOF_estimation, floor(gsPointCount * (9/16)))
    
                currentTime = time;
                tof_current = tof;
        
                deltaV = optimalDVSolver([currentTime, tof_current]);
            end
        end
    end


    %All searches are complimented by a final search for the local minimum 
    %close to the best found solution
    solveDate = 1;
    tf_fuelOptimal = fminsearch(@optimalDVSolver, [dateOptimal, tof_optimal], opt_dv_fminsearch);
    fprintf("Optimization Completed! ---- Finalized %cV Found: <strong>%.0f m/s</strong> \n", 916, deltaResult);

    %contourf(contourMap, 100);

    %plot3(tw_graph(:,1), tw_graph(:,2), tw_graph(:,3),'o','Color','k','MarkerSize',10)
    title(sprintf("DeltaV Map For Different Launch Dates and Time of Flights"));
    xlabel("Time past initial Time [s]");
    ylabel("Time of Flight [s]")
    
    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    
    thrustCurve_date = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_plot, 1)]; 

    deltaV_opt = trapz(theta_vec_plot, abs(fJerkFunction(d_opt, theta_vec_plot, 1)));

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
