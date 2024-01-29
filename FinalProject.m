%% Shape-Based Approach to Low-Thrust Rendezvous Trajectory Design
tic;
clear;
close all;

%% Declaring globals
global pState;
global paramVector paramVector_opt d_opt dateOptimal tof_optimal deltaResult;

%% Central body information

%Gravitational parameter
mju = 3.986004418*10^14; %Earth
%mju = 1.32712440018*10^20; %Sun

%How large is the central body
bodyRadius = 6371 * 1e3; %[km] Earth
%bodyRadius = 696340 * 1e3; %[km] Sun

%Minimum allowed radius
rMin = bodyRadius + 250 * 1e3; %[km] Earth
%rMin = bodyRadius + 1000000 * 1e3; %[km] Sun

%% Orbital parameters
seed = floor(rand() * 1000000);
rng(seed);  % Cool Seeds: Good looking: 761982    Good looking: 13     Some issue to fix: 34158     NuPredict Issue: 438758      DV optim issues: 313

%Save current date and seed to a file
baseDate = datetime('now','Format','dd-MMM-yyyy');
seedFile = fopen("seedList.txt", "a");
fprintf(seedFile, "Seed: %.0f - Date: %s\n", rng().Seed, string(baseDate));

%--First Orbit Parameters--
%Semimajor axis
a_initial = rMin + (0.75 + 10*rand())*rMin;
%Period of the orbit
P1 = 2*pi/sqrt(mju/a_initial^3);
%Time of last perigee pass
Tp1 = rand() * P1;
%Eccentricity
e1 = rand() * 0.8;
%Argument of perigee
omega1 = rand() * 2 * pi;

%--Second Orbit Parameters--
%Semimajor axis
a_final = rMin + (0.75 + 10*rand())*rMin;
%Period of the orbit
P2 = 2*pi/sqrt(mju/a_final^3);
%Time of last perigee pass
Tp2 = rand() * P2;
%Eccentricity
e2 = rand() * 0.8;
%Argument of perigee
omega2 = rand() * 2 * pi;

%% Program settings
%Number of additional rotations around central body
N = 0;

%Are unsolvable positions filled in with maxDepthN options
useMultiorbitFilling = 1;
%What is the maximum depth searched to
maxDepthN = 2;

%What is the average TOF for the orbit
TOF_average = (2*N + 1)*pi*sqrt((a_initial+a_final)^3/(8*mju));

%Limits for TOF in relation to first guess
TofLowMult = 0.2;
TofHighMult = 6;
TofLimLow = TofLowMult * TOF_average;
TofLimHigh = TofHighMult * TOF_average;
%--Adjustment Parameters--
%How much the calculated TOF is increased
TOF_corrMult = 1;
%How much d-coefficients are changed to find new solutions
dAdjustment = 1.025;

%Multiplier to increase the geometric minimum angle
safeTransferAngleMultiplier = 1;

%Maximum allowed radius
rMax = 10*max(a_initial, a_final);

%Spacecraft mass [16U CubeSat]
m = 32; %kg

%Solve TOF, optimize deltaV and/or optimize transfer date
optimizeTOF = 0;
optimizeDV = 0;
optimizeDATE = 0;
optimizeSwarm = 1;

%Accuracies of approximation
intApprox = 100;
plotAccuracy = 500;
fzero_d_accuracy = 10^floor(log10(0.0001 / max(a_initial*(1+e1), a_final*(1+e2))));
fzero_f_accuracy = min(P1, P2) / 1000000;

%Doesn't change, but here not to be a hardcoded value
theta_0 = 0;

%Seconds forwrd from current date
initialTime = 0;

%Impossible time to trigger parameter generation on first run
%if previousTime ~= currentTime
previousTime = -1;

%How far into the future optimal date is searched.
%Testing out lcm to always find interesting values
years1 = ceil(P1/(365*86400));
years2 = ceil(P2/(365*86400));
dateSearchSpan = 3 * max(P1,P2); %0.5 * lcm(years1, years2) * 365 * 86400;

%How many initial points the global search starts with
%Linear for option 1
%Linear for option 2
%Squared for option 3
gsPointCount = 32;

%Which approach to global search is taken
%Option 1: Global search
%Option 2: TOF search with date vector
%Option 3: TOF and date vectors
transferWindowSearchOption = 3;

%How many start TOFs option 2 tries (even 1->2 improves accuracy tremendously)
option2starts = 3;

%How many more points are started with in option 1
option1pointMultiplier = 3;

%Is the transfer window plotted
visualizeTransferWindow = 1;

%Is the Swarm Transfer animation plotted
visualizeSwarmTransfer = 0;

%Is the Swarm Transfer saved as a video
%[visualizeSwarmTransfer] needs to be turned on
saveSwarmTransferVideo = 0;

%Set up the initial deltaV value for plotting
initial_DeltaV = 1e24; %A big number
resultsDeltaVs = [];
%Predefine thrust curves
thrustCurve_tof = [];
thrustCurve_dV = [];
thrustCurve_date = [];

%Configure subplots
subXCount = optimizeTOF + optimizeDV + optimizeDATE;
if subXCount > 0
    subYCount = 5;
    figure(1);
    infoWindow = tiledlayout(subYCount, subXCount*3 + 2, 'Padding', 'tight', 'TileSpacing', 'tight');
    sgtitle(sprintf("Initial Time: %s - Seed: %0.f", secToTime(initialTime, 1, baseDate), rng().Seed));
    
    %Maximize window
    set(gcf,'WindowState','maximized')
end

if (visualizeTransferWindow == 1) && (optimizeDATE == 1)
    % Orbital viewport figure
    figure(3);
    windowSize = get(gcf, 'Position');
    windowSize(3) = windowSize(3) * 0.6;
    set(gcf, 'Position', windowSize)

    viewportWindow = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
    nexttile(viewportWindow);

    xlabel("Distance from central body [m]");
    ylabel("Distance from central body [m]");

    axHandle = gca;    
    hold(axHandle, 'on');
end

%Is the live plotting running
isRunning = 0;

%Create colormap
cbResolution = 100;
cmap = zeros(100, 3);
dv_v = linspace(1, 0, cbResolution);
dv_i = 0.5;

cmap(1,:) = [0.90, 0, 0];

for i = 2:cbResolution
    dv = dv_v(i);
    R_Multiplier = (3.*(dv./dv_i).^2 - 2.*(dv./dv_i).^3);
    G_Multiplier = 1-(3.*((dv-dv_i)./(dv_i)).^2 - 2.*((dv-dv_i)./(dv_i)).^3);
    
    if dv > 2*dv_i
        color = [1 0 0];
    elseif dv > dv_i
        color = [1 sqrt(G_Multiplier) 0];
    elseif dv <= dv_i
        color = [R_Multiplier^2 1 0];
    end
    
   cmap(i,:) = color;
end

%Create plotting data matrix
twMap = zeros(gsPointCount);
twMapInds = [1, 1];

%How many labels does the transfer window have
twLabelCount = 10;

%How thick are the lines in the deltaV bar plot
dvLineWidth = 40;

%% Define some optimization options here for speeeeeeeeed!
opt_tof_fzero = [fzero_d_accuracy, fzero_f_accuracy]; % This is different due to the custom fzero function
opt_tof_fzero_acc = optimset('TolX', 1e-18);%, 'TolFun', 1e-15, 'TolCon', 1e-15, 'TolX', 1e-15);
opt_nu_fzero = optimset('TolFun', 1e-3, 'TolX', 1e-3, 'Display', 'off');
opt_tf_angle = optimset('TolFun', 1e-3, 'Display', 'off');
opt_d_lim_fzero = optimset('TolFun', 1e2, 'TolX', 1e-15, 'Display', 'off');
opt_dv_fminsearch = optimset('TolFun',1e1, 'TolX', min(P1, P2)/1000);
opt_dv_global = optimset('TolFun',1e2, 'TolX', 1e4);
opt_minTheta_fbnd = optimset('TolFun', 1e-2, 'TolX', 1e-2);
opt_d_bounds_fzero = optimset('TolFun', 1e-3);
opt_d_thetaTime_fzero = optimset('TolX', 1e-4);

%% Calculating the orbital parameters
n1 = sqrt(mju/a_initial^3);
p1 = a_initial * (1-e1^2);

n2 = sqrt(mju/a_final^3);
p2 = a_final * (1-e2^2);

%How large are the pixels in the transfer window plot
tfWindowPixelsX = dateSearchSpan / (gsPointCount-1);
tfWindowPixelsY = (TofLimHigh - TofLimLow) / (gsPointCount-1);

%Propagating the two orbits
nu = linspace(0, 2*pi, plotAccuracy);
orbit1 = [cos(nu+omega1) * p1 ./ (1+e1*cos(nu)); sin(nu+omega1) * p1 ./ (1+e1*cos(nu))];
orbit2 = [cos(nu+omega2) * p2 ./ (1+e2*cos(nu)); sin(nu+omega2) * p2 ./ (1+e2*cos(nu))];

%Update viewport limits based on orbits
allX = [orbit1(1,:), orbit2(1,:)];
allY = [orbit1(2,:), orbit2(2,:)];
extraDist = (a_final + a_initial) / 4;
portSideLen = max([allX allY]) - min([allX allY]) + 2*extraDist;
portCenter = [(min(allX) + max(allX)) / 2, (min(allY) + max(allY)) / 2];
%viewportLimits = [min([allX allY]) - extraDist, max([allX allY]) + extraDist, min([allX allY]) - extraDist, max([allX allY]) + extraDist];
viewportLimits = [portCenter(1)-portSideLen / 2, portCenter(1)+portSideLen / 2, portCenter(2)-portSideLen / 2, portCenter(2)+portSideLen / 2];    

%% Create and update global setting and state classes

pSettings = programSettings;
pState = programState;
paramVector = paramClass;
paramVector.mju = mju;
paramVector.theta1 = 0;
paramVector.theta2 = 0;

pSettings.mju = mju;
pSettings.rMin = rMin;
pSettings.rMax = rMax;
pSettings.a_initial = a_initial;
pSettings.a_final = a_final;
pSettings.omega1 = omega1;
pSettings.omega2 = omega2;
pSettings.e1 = e1;
pSettings.e2 = e2;
pSettings.Tp1 = Tp1;
pSettings.Tp2 = Tp2;
pSettings.n1 = n1;
pSettings.n2 = n2;
pSettings.P1 = P1;
pSettings.P2 = P2;
pSettings.p1 = p1;
pSettings.p2 = p2;
pSettings.intApprox = intApprox;
pSettings.approxSpace = linspace(0, 1, intApprox);
pSettings.plotAccuracy = plotAccuracy;
pSettings.theta_0 = theta_0;
pSettings.safeTransferAngleMultiplier = safeTransferAngleMultiplier;
pSettings.TOF_corrMult = TOF_corrMult;
pSettings.TOF_average = TOF_average;
pSettings.dAdjustment = dAdjustment;
pSettings.opt_nu_fzero = opt_nu_fzero;
pSettings.opt_tf_angle = opt_tf_angle;
pSettings.opt_tof_fzero = opt_tof_fzero;
pSettings.opt_d_lim_fzero = opt_d_lim_fzero;
pSettings.opt_minTheta_fbnd = opt_minTheta_fbnd;
pSettings.opt_d_bounds_fzero = opt_d_bounds_fzero;
pSettings.opt_dv_fminsearch = opt_dv_fminsearch;
pSettings.opt_tof_fzero_acc = opt_tof_fzero_acc;
pSettings.opt_dv_global = opt_dv_global;
pSettings.solveDate = 0; %Default as 0
pSettings.plotTransferWindow = 0; %Default as 0
pSettings.useMultiorbitFilling = useMultiorbitFilling;
pSettings.tfWindowPixelsX = tfWindowPixelsX;
pSettings.tfWindowPixelsY = tfWindowPixelsY;
pSettings.maxDepthN = maxDepthN;
pSettings.gsPointCount = gsPointCount;
pSettings.twLabelCount = twLabelCount;
pSettings.transferWindowSearchOption = transferWindowSearchOption;
pSettings.option2starts = option2starts;
pSettings.option1pointMultiplier = option1pointMultiplier;
pSettings.cmap = cmap;
pSettings.bodyRadius = bodyRadius;
pSettings.orbit1 = orbit1;
pSettings.orbit2 = orbit2;
pSettings.viewportLimits = viewportLimits;
pSettings.baseDate = baseDate;

pState.currentTime = initialTime;
pState.tof_current = 0; %Default as 0
pState.initial_tof = 0; %Default as 0
pState.previousTime = previousTime;
pState.N = N;
pState.initial_DeltaV = initial_DeltaV;
pState.successfulOrbits = 0; %Default as 0
pState.testedOrbits = 0; %Default as 0
pState.isRunning = isRunning;
pState.axRelative = [0, 0, 0, 0]; %Default as 0
pState.twMap = twMap;
pState.twMapInds = twMapInds;

%% Calculate all the rest of orbital parameters for the coming simulation
[resultVector, paramVector] = updateParameters(1, pSettings);

theta_f = paramVector.theta_f;
r1 = paramVector.r1;
r2 = paramVector.r2;
theta1 = paramVector.theta1;
theta2 = paramVector.theta2;
nu2_i = paramVector.nu2_i;
r2_i = paramVector.r2_i;

d_minimum = resultVector(1);
d_maximum = resultVector(2);
TOF_estimation = pState.initial_tof;

%% TOF optimization + result plotting
if optimizeTOF == 1 
    solutionFound = 0;
    startN = pState.N;

    %Plot results of TOF solved trajectory
    nexttile(infoWindow, [3, 3]);
    hold on;
    while (solutionFound == 0)
        try
            %TOF solution
            theta_vec_plot = linspace(theta_0, theta_f, plotAccuracy);
            theta_vec_super = fGetThetaSuper(theta_vec_plot);
            dT = theta_vec_super(1,2) - theta_vec_super(1,1);

            tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, pState.tof_current, theta_vec_super);
            d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], opt_tof_fzero_acc);
            solutionFound = 1;
            break;
        catch

            pState.N = pState.N + 1;
            if pState.N > pSettings.maxDepthN
                break
            end
            [resultVector, paramVector] = updateParameters(1, pSettings);

        end
    end

    pState.N = startN;

    
    xlabel("Distance from central body [m]");
    ylabel("Distance from central body [m]");

    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
        
    rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
        
    if solutionFound == 1

        r1 = paramVector.r1;
        r2 = paramVector.r2;
        theta1 = paramVector.theta1;
        theta2 = paramVector.theta2;
        nu2_i = paramVector.nu2_i;
        r2_i = paramVector.r2_i;
      
        d_minimum = resultVector(1);
        d_maximum = resultVector(2);
        TOF_estimation = pState.initial_tof;

        thrustCurve_tof = [theta_vec_plot; fThrustFunction(d_solution, theta_vec_super, paramVector)];    
    
        deltaV_tof = fJerkFunction(d_solution, theta_vec_super, dT, paramVector);
        pState.initial_DeltaV = deltaV_tof;
        resultsDeltaVs(end+1) = deltaV_tof;
        
        time_t = fTimeFunction(d_solution, theta_vec_super, dT, paramVector);
        
        plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
        plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')   
        
        x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_super, paramVector);
        y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_solution, theta_vec_super, paramVector);
        plot(x, y, "Color", [0.2 0.7 0.2]);

        title(sprintf("Estimated TOF solution\nTarget TOF: %s\nAchieved TOF: %s\nRequired deltaV: %.0f m/s", secToTime(TOF_estimation, 0, []), secToTime(time_t, 0, []), deltaV_tof));
   
    else
        fprintf("Inital Time of Flight guess not achievable\n")

        %Reset revolutions
        pState.N = startN;
        [resultVector, paramVector] = updateParameters(1, pSettings);

        theta_f = paramVector.theta_f;
        theta1 = paramVector.theta1;
        theta2 = paramVector.theta2;
        nu2_i = paramVector.nu2_i;
        r2 = paramVector.r2;
        r2_i = paramVector.r2_i;
        
        d_minimum = resultVector(1);
        d_maximum = resultVector(2);
        TOF_estimation = pState.initial_tof;

        theta_vec_plot = linspace(0, theta_f, plotAccuracy);
        theta_vec_super = fGetThetaSuper(theta_vec_plot);
        dT = theta_vec_super(1,2) - theta_vec_super(1,1);

        plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
        plot(cos(omega2 + nu2_i) * r2_i, sin(omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')   
        
        x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_minimum, theta_vec_super, paramVector);
        y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_minimum, theta_vec_super, paramVector);
        time_max = fTimeFunction(d_minimum, theta_vec_super, dT, paramVector);
        plot(x, y, "Color", [0.5 0.9 0.5]);

        x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_maximum, theta_vec_super, paramVector);
        y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_maximum, theta_vec_super, paramVector);
        time_min = fTimeFunction(d_maximum, theta_vec_super, dT, paramVector);
        plot(x, y, "Color", [0.5 0.9 0.5]);

        title(sprintf("Inital Time of Flight guess not achievable\nTarget TOF: %s\nMinimum Achieved TOF: %s\nMaximum Achieved TOF: %s", secToTime(TOF_estimation, 0, []), secToTime(time_min, 0, []), secToTime(time_max, 0, [])));
        
        thrustCurve_tof = [0;0];

    end

    %legend("Initial orbit", "Target orbit", "", "", "", "Transfer Orbit");

    xlim(viewportLimits(1:2));
    ylim(viewportLimits(3:4));
end
%% Delta V optimization for fixed starting point + result plotting
if optimizeDV == 1
    pSettings.solveDate = 0;
    pSettings.plotTransferWindow = 0;

    %Initialize best values
    deltaResult = Inf;

    dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
    
    for tof = [TofLimHigh, TOF_estimation, TofLimHigh]
        tf_fuelOptimal = fminsearch(dvHandle, tof, opt_dv_fminsearch);
    end

    %Update local variables to the optimal solution
    theta_f_opt = paramVector_opt.theta_f;
    r1_opt = paramVector_opt.r1;
    r2_opt = paramVector_opt.r2;
    theta1_opt = paramVector_opt.theta1;
    theta2_opt = paramVector_opt.theta2;
    nu2_i_opt = paramVector_opt.nu2_i;
    r2_i_opt = paramVector_opt.r2_i;

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    theta_vec_super = fGetThetaSuper(theta_vec_plot);
    dT = theta_vec_super(1,2) - theta_vec_super(1,1);
    
    thrustCurve_dV = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_super, paramVector_opt)]; 
    
    deltaV_opt = fJerkFunction(d_opt, theta_vec_super, dT, paramVector_opt);
    pState.initial_DeltaV = deltaV_opt;
    resultsDeltaVs(end+1) = deltaV_opt;

    time_t = fTimeFunction(d_opt, theta_vec_super, dT, paramVector_opt);
    
    nexttile(infoWindow, [3, 3]);
    hold on;
    
    xlabel("Distance from central body [m]");
    ylabel("Distance from central body [m]");
    
    plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")

    plot(cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
        
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_super, paramVector_opt);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_super, paramVector_opt);
    plot(x, y, "Color", [0.2 0.7 0.2]);
        
    title(sprintf("DeltaV optimized TOF\nTransfer date: %s\nUsed TOF: %s\nRequired deltaV: %.0f m/s", secToTime(initialTime, 1, baseDate), secToTime(time_t, 0, []), deltaV_opt));

    %legend("Initial orbit", "Target orbit", "", "", "", "Transfer Orbit");

    xlim(viewportLimits(1:2));
    ylim(viewportLimits(3:4));
end
%% Optimize the transfer date for deltaV + result plotting
if optimizeDATE == 1
    pSettings.plotTransferWindow = visualizeTransferWindow;

    %Initialize best values
    deltaResult = Inf;
    
    %Reset current Time
    pState.tof_current = TOF_estimation;

    %Befor changing figure
    nexttile(infoWindow, [3 3]);
    xlabel("Distance from central body [m]");
    ylabel("Distance from central body [m]");
    orbitAx = gca;

    if visualizeTransferWindow == 1
        twFig = figure(2);
        transferWindow = tiledlayout(1, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
        nexttile;

        % Set up the zoom callback function
        set(zoom, 'ActionPostCallback', @(src, event) zoomCallback(src, event, pSettings));
    
        hold on;
        %Maximize window
        set(gcf,'WindowState','maximized')
        title(sprintf("DeltaV Map For Different Launch Dates and Time of Flights - Seed: %.0f", rng().Seed));
        xlabel("Launch date");
        ylabel("Time of Flight");

        colormap(cmap)
        cb = colorbar;
        clim([0, 1])
        cb.Ticks = [0, 0.5, 1];
        cb.TickLabels = [{sprintf('More %cV', 916)}, {sprintf('Same %cV', 916)}, {sprintf('Less %cV', 916)}];
        
        xlim([initialTime, initialTime + dateSearchSpan])
        ylim([TofLimLow, TofLimHigh])       

        tof_vec = linspace(TofLimLow, TofLimHigh, pSettings.twLabelCount);
        ylims = ylim;
        tof_coords = linspace(ylims(1), ylims(2), pSettings.twLabelCount);

        date_vec = linspace(initialTime, initialTime + dateSearchSpan, pSettings.twLabelCount);
        xlims = xlim;
        date_coords = linspace(xlims(1), xlims(2), pSettings.twLabelCount);

        xticks(date_coords);
        xticklabels(string(secToTime(date_vec, 1, baseDate)));

        yticks(tof_coords);
        tickTexts = cell(1, pSettings.twLabelCount);
        for i = 1:pSettings.twLabelCount
            tickTexts(i) =  {sprintf("%.2f h", tof_vec(i) / (3600))};
        end
        yticklabels(tickTexts);
    end

    pState.successfulOrbits = 0;
    pState.testedOrbits = 0;
    %-- Solve transfer window using Global Search --
    if transferWindowSearchOption == 1
        pSettings.solveDate = 1;
        gs = MultiStart;
        numberOfStartPoints = gsPointCount * option1pointMultiplier;
        dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
        problem = createOptimProblem('fmincon','x0',[initialTime, TOF_average],...
        'objective',dvHandle,'lb',[initialTime, TofLimLow], ...
        'ub',[initialTime + dateSearchSpan, TofLimHigh], 'options', opt_dv_global);
        run(gs, problem, numberOfStartPoints);
    end    
    
    %-- Solve transfer window using set dates and optimize TOF --
    if transferWindowSearchOption == 2
        for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong> \n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof_start = linspace(TofLimLow + (TofLimHigh - TofLimLow)/option2starts, TofLimHigh, option2starts)
                pSettings.solveDate = 0;
                pState.currentTime = time;
                dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
                fminsearch(dvHandle, tof_start, opt_dv_fminsearch);
            end
        end
    end

    %-- Solve transfer window using a set grid of dates and TOFs --
    if transferWindowSearchOption == 3
        pSettings.solveDate = 1;
        for time = linspace(initialTime, initialTime + dateSearchSpan, gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong>\n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof = linspace(TofLimLow, TofLimHigh, gsPointCount)
                pState.currentTime = time;
                pState.tof_current = tof;
                optimalDVSolver([pState.currentTime, pState.tof_current], pSettings, pState.twMapInds);
                pState.twMapInds(2) = pState.twMapInds(2) + 1;
            end
            pState.twMapInds(1) = pState.twMapInds(1) + 1;
            pState.twMapInds(2) = 1;
        end
        %Plot the image
        imagesc([initialTime, initialTime + dateSearchSpan], [TofLimLow, TofLimHigh], pState.twMap);
        %colormap(cmap);
    end
    
    %All searches are complimented by a final search for the local minimum 
    %close to the best found solution
    pSettings.solveDate = 1;
    pSettings.plotTransferWindow = 0;
    dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
    tf_fuelOptimal = fminsearch(dvHandle, [dateOptimal, tof_optimal], opt_dv_fminsearch);
    fprintf("Optimization Completed! ---- Finalized %cV Found: <strong>%.0f m/s</strong> \n", 916, deltaResult);

    fprintf("Tested %.0f transfer orbits out of which %1.f %% (%.0f) were achievable\n", pState.testedOrbits, 100 * (pState.successfulOrbits / pState.testedOrbits), pState.successfulOrbits);

    if (visualizeTransferWindow == 1) && (transferWindowSearchOption ~= 3)
        rectangle('Position',[dateOptimal-0.5*pSettings.tfWindowPixelsX, tof_optimal-0.5*pSettings.tfWindowPixelsY, ...
                   pSettings.tfWindowPixelsX, pSettings.tfWindowPixelsY], 'LineWidth', 1, 'EdgeColor', 'black');
    end

    %Update local variables to the optimal solution
    theta_f_opt = paramVector_opt.theta_f;
    r1_opt = paramVector_opt.r1;
    r2_opt = paramVector_opt.r2;
    theta1_opt = paramVector_opt.theta1;
    theta2_opt = paramVector_opt.theta2;
    nu2_i_opt = paramVector_opt.nu2_i;
    r2_i_opt = paramVector_opt.r2_i;

    theta_vec_plot = linspace(theta_0, theta_f_opt, plotAccuracy);
    theta_vec_super = fGetThetaSuper(theta_vec_plot);
    dT = theta_vec_super(1,2) - theta_vec_super(1,1);

    thrustCurve_date = [theta_vec_plot; fThrustFunction(d_opt, theta_vec_super, paramVector_opt)]; 

    deltaV_date = fJerkFunction(d_opt, theta_vec_super, dT, paramVector_opt);
    resultsDeltaVs(end+1) = deltaV_date;

    time_t = fTimeFunction(d_opt, theta_vec_super, dT, paramVector_opt);

    hold(orbitAx, 'on') 
    
    plot(orbitAx, orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
    plot(orbitAx, orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
    
    plot(orbitAx, cos(theta1_opt) * r1_opt, sin(theta1_opt) * r1_opt,'or', 'MarkerSize',5,'MarkerFaceColor','g')
    plot(orbitAx, cos(theta2_opt) * r2_opt, sin(theta2_opt) * r2_opt,'or', 'MarkerSize',5,'MarkerFaceColor','r')
    plot(orbitAx, cos(omega2 + nu2_i_opt) * r2_i_opt, sin(omega2 + nu2_i_opt) * r2_i_opt,'or', 'MarkerSize',5,'MarkerFaceColor','k')
        
    rectangle(orbitAx, 'Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
    
    x = cos(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_super, paramVector_opt);
    y = sin(theta_vec_plot+theta1_opt) .* fRadiusFunction(d_opt, theta_vec_super, paramVector_opt);
   
    plot(orbitAx, x, y, "Color", [0.2 0.7 0.2]);

    title(orbitAx, sprintf("Full transfer window solution\nTransfer date: %s\nUsed TOF: %s\nRequired deltaV: %.0f m/s",secToTime(dateOptimal, 1, baseDate), secToTime(time_t, 0, []), deltaV_date));
    
    %legend(orbitAx, "Initial orbit", "Target orbit", "", "", "", "Transfer Orbit");
    
    orbitAx.XLim = viewportLimits(1:2);
    orbitAx.YLim = viewportLimits(3:4);
end

%% Find solutions for a swarm of satellites
if optimizeSwarm
    %global resultVector_opt

    leadNuPointCountPerOrbit = 25;
    datePointCount = 25;

    extraOrbitChecks = 0;

    %Very primitive logic in place.
    extraOrbitChecks = extraOrbitChecks + floor((P1/P2)^(1/3));

    satelliteCount = 10;

    deployTime = initialTime;
    targetTime = 0;

    %What are the true anomalies of the objects after deployment?
    deployNuVector = linspace(0, 0.1, satelliteCount) .* (2 * pi);

    %Needs to be generated from a time vector to get true equal spacing on
    %elliptical orbits
    perigeeTimeVectorFinal = linspace(0, (1-1/satelliteCount) * pSettings.P2, satelliteCount);

    targetNuVector = zeros(1, satelliteCount);
    for ind = 1:satelliteCount
        T_nu = perigeeTimeVectorFinal(ind);
        if T_nu ~= 0        
            n = pSettings.n2;
            e = pSettings.e2;
            nu_guess = 0;
    
            nu2 = nuFromTime(T_nu, n, e, nu_guess);
        else
            nu2 = 0;
        end
        targetNuVector(ind) = nu2;
    end
   
    nuCount = size(deployNuVector);
    nuCount = nuCount(2);

    pSettings.solveDate = 1;
    pSettings.plotTransferWindow = 0;
    pSettings.useMultiorbitFilling = 0;

    perigeeTimeVectorInitial = zeros(1,nuCount);

    TpVectorInd = 1;

    for nuStart = deployNuVector
        if nuStart <= pi
            E = 2*atan(tan(nuStart/2)/sqrt((1+pSettings.e1)/(1-pSettings.e1)));
        else
            E = 2*pi + 2*atan(tan(nuStart/2)/sqrt((1+pSettings.e1)/(1-pSettings.e1)));
        end
        perigeeTimeVectorInitial(TpVectorInd) = (E-pSettings.e1*sin(E)) / pSettings.n1;
        TpVectorInd = TpVectorInd + 1;
    end   

    tofMatrix = zeros(nuCount, datePointCount, leadNuPointCountPerOrbit * (extraOrbitChecks+1));
    nuInd = 1;
    tofMatDateInd = 1;
    leadNuIndex = 1;

    dateVector = linspace(initialTime, initialTime + dateSearchSpan, datePointCount);
   
    %There needs to be logic to determine how much lead time is given based
    %on relative orbit sizes. 
    leadTimeVector = linspace(-pSettings.P2, pSettings.P2, leadNuPointCountPerOrbit);
        
    for targetPerigeeTime = perigeeTimeVectorFinal
        for date = dateVector
            for leadOffset = leadTimeVector
                leadTime = leadOffset - perigeeTimeVectorFinal(nuInd);
                requiredTOF = mod(targetTime - date, pSettings.P2) + leadTime;
                tofVector = requiredTOF + (0:extraOrbitChecks) .* pSettings.P2;
                tofMatrix(nuInd, tofMatDateInd, (leadNuIndex : leadNuIndex + extraOrbitChecks)) = tofVector;
                leadNuIndex = leadNuIndex + (extraOrbitChecks+1);
            end
            tofMatDateInd = tofMatDateInd + 1;
            leadNuIndex = 1;
        end   
        nuInd = nuInd + 1;
        tofMatDateInd = 1;
    end
 
    pState.currentTime = initialTime;

    swarmParams = swarmParamClass;

    swarmParams.leadNuPointCountPerOrbit = leadNuPointCountPerOrbit;
    swarmParams.datePointCount = datePointCount;
    swarmParams.extraOrbitChecks = extraOrbitChecks;
    swarmParams.nuCount = nuCount;
    swarmParams.dateVector = dateVector;
    swarmParams.tofMatrix = tofMatrix;
    swarmParams.deployTime = deployTime;
    swarmParams.targetTime = targetTime;

    swarmParams.perigeeTimeVectorInitial = perigeeTimeVectorInitial;
    swarmParams.perigeeTimeVectorFinal = perigeeTimeVectorFinal;
    
    %Find best solution for each combination orbit in the swarm
    [resultMatrix, dateResultMatrix, tofResultMatrix, paramResultMatrix, efgResultMatrix] = optimalSwarmDVSolver(swarmParams, pSettings,[]);
    
    %Find the optimal combination of the orbits using the Hungarian /
    %Munkres Algorithm. Not my implementation. Credit to Yi Cao at Cranfield University
    [assignment, cost] = fMunkresAlgorithm(resultMatrix);
    
    sortedInds = sortrows([assignment; (1:nuCount)], 1);
    assignment = sortedInds(2,:);

    requiredDeltaVs = zeros(nuCount, 1);
    requiredDates = zeros(nuCount, 1);
    requiredTOFs = zeros(nuCount, 1);
    
    for resultInd = 1:nuCount
        requiredDeltaVs(resultInd) = resultMatrix(resultInd, sortedInds(resultInd));
        requiredDates(resultInd) = dateResultMatrix(resultInd, sortedInds(resultInd));
        requiredTOFs(resultInd) = tofResultMatrix(resultInd, sortedInds(resultInd));
    end

    totalDV = sum(requiredDeltaVs);
    
    fprintf("Average Required %cV: <strong>%.0f m/s</strong>, First Launch: %s, Last Launch: %s\n", 916, totalDV/satelliteCount, secToTime(min(requiredDates), 1, pSettings.baseDate), secToTime(max(requiredDates), 1, pSettings.baseDate))

    for deployInd = 1:nuCount
        targetInd = assignment(deployInd);
        fprintf("Object %.0f:\tRequired: <strong>%cV %.0f m/s</strong>, Launch date: %s, TOF: %s\n", deployInd, 916, resultMatrix(deployInd, targetInd), secToTime(dateResultMatrix(deployInd, targetInd), 1, pSettings.baseDate), secToTime(tofResultMatrix(deployInd, targetInd), 0, []));        
    end

    if visualizeSwarmTransfer == 1
    
        % Define the desired axis size in pixels
        axisWidthPixels = 540;
        axisHeightPixels = 540;
        swarmPlot = figure(4);
        set(swarmPlot, 'Units', 'pixels');
        set(swarmPlot, 'Position', [40, 40, axisWidthPixels + 40, axisHeightPixels + 40]);
    
        % Get the figure size in pixels
        figurePosition = get(gcf, 'Position');
        figureWidthPixels = figurePosition(3);
        figureHeightPixels = figurePosition(4);
        
        % Calculate the normalized axis size
        axisWidthNormalized = axisWidthPixels / figureWidthPixels;
        axisHeightNormalized = axisHeightPixels / figureHeightPixels;
        
        % Calculate the normalized axis position (assuming it starts at the lower-left corner)
        axisLeftNormalized = 0.5 * (figureWidthPixels - axisWidthPixels) / figureWidthPixels;
        axisBottomNormalized = 0.5 * (figureHeightPixels - axisHeightPixels) / figureHeightPixels;
        
        % Set the axis size in pixels
        axes('Position', [axisLeftNormalized, axisBottomNormalized, axisWidthNormalized, axisHeightNormalized]);
    
        hold on;
    
        xlabel("Distance from central body [m]");
        ylabel("Distance from central body [m]");
    
        plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
        plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2)
    
        axis manual
    
        if saveSwarmTransferVideo == 1
            % Create a VideoWriter object
            writerObj = VideoWriter('animated_video.mp4', 'MPEG-4');
            writerObj.FrameRate = 60;  % Set the frame rate (frames per second)
            
            % Open the video file for writing
            open(writerObj);
        end
    
        interParam = paramClass;
        efg_Mat_1 = zeros(3,3);
    
        swarmTimeSpan = max(requiredTOFs) + max(requiredDates - deployTime);
        firstDeployTime = min(requiredDates);
    
        for time = linspace(firstDeployTime - pSettings.P1, (firstDeployTime + swarmTimeSpan) + pSettings.P2, 2*pSettings.plotAccuracy)
            hold off
            plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
            hold on
            plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
            rectangle('Position',[-bodyRadius, -bodyRadius, 2*bodyRadius, 2*bodyRadius],'Curvature',[1 1], 'FaceColor',"yellow")
    
            for deployInd = 1:nuCount
                targetInd = assignment(deployInd);
    
                if time < (dateResultMatrix(deployInd, targetInd))% + tofResultMatrix(deployInd, targetInd))
                    pState.tof_current = 0;
                    pState.currentTime = time;
    
                    Tp1 = perigeeTimeVectorInitial(deployInd);
                    Tp2 = perigeeTimeVectorFinal(targetInd);
                    [~, paramVector, ~] = updateSwarmParameters(Tp1, Tp2, deployTime, targetTime, 1, [], [], pSettings);
    
                    plot(cos(paramVector.theta1) * paramVector.r1, sin(paramVector.theta1) * paramVector.r1,'or', 'MarkerSize',5,'MarkerFaceColor','g');
    
                elseif time > (dateResultMatrix(deployInd, targetInd) + tofResultMatrix(deployInd, targetInd))
                    pState.tof_current = 0;
                    pState.currentTime = time;
    
                    Tp1 = perigeeTimeVectorInitial(deployInd);
                    Tp2 = perigeeTimeVectorFinal(targetInd);
                    [~, paramVector, ~] = updateSwarmParameters(Tp1, Tp2, deployTime, targetTime, 1, [], [], pSettings);
    
                    plot(cos(paramVector.theta2) * paramVector.r2, sin(paramVector.theta2) * paramVector.r2,'or', 'MarkerSize',5,'MarkerFaceColor','r');
                
                else
                    paramResultVector = paramResultMatrix(deployInd, targetInd, :);
    
                    efgVector = reshape(efgResultMatrix(deployInd, targetInd, :), [1,9]);
                    efg_Mat_1_result = [efgVector(1:3); efgVector(4:6); efgVector(7:9)];
                    
                    interParam.mju =  paramResultVector(1);
                    interParam.gamma1 =  paramResultVector(2);
                    interParam.gamma2 =  paramResultVector(3);
                    interParam.theta_f =  paramResultVector(4);
                    interParam.theta1_dot =  paramResultVector(5);
                    interParam.theta2_dot =  paramResultVector(6);
                    interParam.r1 =  paramResultVector(7);
                    interParam.r2 =  paramResultVector(8);
                    interParam.theta1 =  paramResultVector(9);
                    interParam.theta2 = paramResultVector(10);
                    interParam.a = paramResultVector(11);
                    interParam.b = paramResultVector(12);
                    interParam.c = paramResultVector(13);
                    interParam.efg_Mat_1 = efg_Mat_1_result;
    
                    d_result = paramResultVector(14);
    
                    transferTimeRatio = time - dateResultMatrix(deployInd, targetInd);
                    normalRatio = transferTimeRatio / tofResultMatrix(deployInd, targetInd);
    
                    fThetaTimeHandle = @(theta) fThetaFromTime(theta, d_result, interParam, transferTimeRatio);
                    
                    %optimset not specific to this function
                    currentTheta = fzero(fThetaTimeHandle, [-0.01, 1.01 * interParam.theta_f], opt_d_thetaTime_fzero);
                    superRatio = fGetThetaSuper(currentTheta);
                    
                    plotPointCount = ceil((1-normalRatio) * plotAccuracy);
                    if plotPointCount > 1
                        theta_vec_plot = linspace(currentTheta, interParam.theta_f, plotPointCount);
                        theta_vec_super = fGetThetaSuper(theta_vec_plot);
                        dT = theta_vec_super(1,2) - theta_vec_super(1,1);
                
                        x = cos(theta_vec_plot+interParam.theta1) .* fRadiusFunction(d_result, theta_vec_super, interParam);
                        y = sin(theta_vec_plot+interParam.theta1) .* fRadiusFunction(d_result, theta_vec_super, interParam);
                        plot(x, y, "Color", [0.2 0.7 0.2]);
                    end
    
                    x = cos(currentTheta + interParam.theta1) .* fRadiusFunction(d_result, superRatio, interParam);
                    y = sin(currentTheta + interParam.theta1) .* fRadiusFunction(d_result, superRatio, interParam);
                    plot(x, y,'or', 'MarkerSize',5,'MarkerFaceColor','y');
    
    
    
                end
            end
            xlim(viewportLimits(1:2));
            ylim(viewportLimits(3:4));
    
            drawnow;
    
            if saveSwarmTransferVideo == 1
                % Capture the frame for the video
                frame = getframe(gcf);
                
                % Write the frame to the video file
                writeVideo(writerObj, frame);
            end
    
        end
    
        if saveSwarmTransferVideo == 1
            close(writerObj);
        end
    end
end



%% Plotting thrust curves and deltaV results
if subXCount > 0
    nexttile(infoWindow, [subYCount, 2]);
    hold on
    ylabel("deltaV [km/s]");
    xticks([]);
    xticklabels([]);
    deltaInd = 1;
    successValues = size(resultsDeltaVs);
    lineXpos = linspace(0, 1, successValues(2));
    newXLim = [lineXpos(1)-0.25, lineXpos(end)+0.25];
    xlim(newXLim);
    
    if (optimizeTOF == 1) && (solutionFound == 1)
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth*1.1, 'Color', 'black');
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth, 'Color', [0.9, 0.8, 0.0]);
        dvtx1 = text(lineXpos(deltaInd), resultsDeltaVs(deltaInd) / 1100, 'Estimated TOF', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90);
        set(dvtx1, 'FontSize', 12, 'FontWeight', 'bold')
    
        deltaInd = deltaInd+1;
    end
    
    if optimizeDV == 1 
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth*1.1, 'Color', 'black');
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth, 'Color', [0.9, 0.0, 0.0]);
        dvtx2 = text(lineXpos(deltaInd), resultsDeltaVs(deltaInd) / 1100, 'Free TOF', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90);
        set(dvtx2, 'FontSize', 12, 'FontWeight', 'bold')
    
        deltaInd = deltaInd+1;
    end
    
    if optimizeDATE == 1
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth*1.1, 'Color', 'black');
        plot([lineXpos(deltaInd), lineXpos(deltaInd)],[0, resultsDeltaVs(deltaInd) / 1000], 'LineWidth', dvLineWidth, 'Color', [0.2, 0.2, 0.9]);
        dvtx3 = text(lineXpos(deltaInd), resultsDeltaVs(deltaInd) / 1100, 'Transfer window', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Rotation', 90);
        set(dvtx3, 'FontSize', 12, 'FontWeight', 'bold')
    end
    
    title(sprintf("Comparision of total deltaV values"));
    
    
    nexttile(infoWindow, [2 subXCount*3])
    hold on;
    xlabel("theta");
    ylabel("thurst [mN]");
    maxThetaValue =  max([thrustCurve_tof(1, end), thrustCurve_dV(1, end), thrustCurve_date(1, end)]);
    plot([0, maxThetaValue], [0, 0], 'Color', 'black', 'LineStyle',':')
    xlim([0, maxThetaValue+0.5])
    
    if (optimizeTOF == 1) && (solutionFound == 1)
        % Plot results of tof optimized trajectory
        plot(thrustCurve_tof(1, :), 1000*m*thrustCurve_tof(2, :), 'LineWidth', 2, 'Color', [0.9, 0.8, 0.0]);
        tctx1 = text(thrustCurve_tof(1, end), 1000*m*thrustCurve_tof(2, end), 'Estimated TOF', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        set(tctx1, 'FontSize', 12, 'FontWeight', 'bold')
    end
    if optimizeDV == 1
        % Plot results of dV optimized trajectory
        plot(thrustCurve_dV(1, :), 1000*m*thrustCurve_dV(2, :), 'LineWidth', 2, 'Color', [0.9, 0.0, 0.0]);
        tctx2 = text(thrustCurve_dV(1, end), 1000*m*thrustCurve_dV(2, end), 'Free TOF', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        set(tctx2, 'FontSize', 12, 'FontWeight', 'bold')
    end
    if optimizeDATE == 1
        % Plot results of date optimized trajectory
        plot(thrustCurve_date(1, :), 1000*m*thrustCurve_date(2, :), 'LineWidth', 2, 'Color', [0.2, 0.2, 0.9]);
        tctx3 = text(thrustCurve_date(1, end), 1000*m*thrustCurve_date(2, end), 'Transfer Window', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        set(tctx3, 'FontSize', 12, 'FontWeight', 'bold')
    end
    
    title(sprintf("Comparision of thrust curves for all solutions"));
end

%% Set up figures and handles
if (visualizeTransferWindow == 1) && (optimizeDATE == 1)    
    % Set up the hover callback function
    set(twFig, 'WindowButtonMotionFcn', @(src, event) hoverCallback(src, event, pSettings, axHandle));
end

toc;

%% START OF FUNCTIONS
%% Second to time conversion
function [timeString] = secToTime(secs, outputDate, baseDate)
    if outputDate == 1
        timeString = baseDate + seconds(secs);

    else
        
        timeStringLimit = 3;
        timeStringCount = 0;
        timeString = "";

        years = floor(secs / (365.25 * 86400));
        secs = secs - years * (365.25 * 86400);
        if (years > 0) && (timeStringCount < timeStringLimit)
            timeString = append(timeString, sprintf("%.0f y ",years));
            timeStringCount = timeStringCount + 1; 
        end

        days = floor(secs / 86400);
        secs = secs - days * 86400;
        if (days > 0) && (timeStringCount < timeStringLimit)
            timeString = append(timeString, sprintf("%.0f d ",days));
            timeStringCount = timeStringCount + 1; 
        end

        hours = floor(secs / 3600);
        secs = secs - hours * 3600;
        if (hours > 0) && (timeStringCount < timeStringLimit)
            timeString = append(timeString, sprintf("%.0f h ",hours));
            timeStringCount = timeStringCount + 1; 
        end

        minutes = floor(secs / 60);
        if (minutes > 0) && (timeStringCount < timeStringLimit)
            timeString = append(timeString, sprintf("%.0f m ",minutes));
            timeStringCount = timeStringCount + 1; 
        end

        secs = secs - minutes * 60;
        if (secs > 0) && (timeStringCount < timeStringLimit)
            timeString = append(timeString, sprintf("%.0f s ",secs));
            timeStringCount = timeStringCount + 1; 
        end

    end

end

%% Create the supervector
function [theta_super] = fGetThetaSuper(theta_vec)
%Pre power the long vector instead of doing it many times
    theta_vec1 = theta_vec;
    theta_vec2 = theta_vec1.*theta_vec;
    theta_vec3 = theta_vec2.*theta_vec;
    theta_vec4 = theta_vec3.*theta_vec;
    theta_vec5 = theta_vec4.*theta_vec;
    theta_vec6 = theta_vec5.*theta_vec;

    theta_super = [theta_vec1; theta_vec2; theta_vec3; theta_vec4; theta_vec5; theta_vec6];
end

%% Zoom callback function
function zoomCallback(~, eventData, pSettings)
    tic;
    global deltaResult;
    global pState;
    pState.isRunning = 1;
 
    %Update transfer window bounds
    initialTime = eventData.Axes.XLim(1);
    dateSearchSpan = eventData.Axes.XLim(2) - initialTime;
    TofLimLow = eventData.Axes.YLim(1);
    TofLimHigh = eventData.Axes.YLim(2);
    %How large are the pixels in the transfer window plot
    tfWindowPixelsX = dateSearchSpan / (pSettings.gsPointCount-1);
    tfWindowPixelsY = (TofLimHigh - TofLimLow) / (pSettings.gsPointCount-1);

    pSettings.tfWindowPixelsX = tfWindowPixelsX;
    pSettings.tfWindowPixelsY = tfWindowPixelsY;
    
    pState.currentTime = initialTime;

    pSettings.plotTransferWindow = 1;

    %Initialize best values
    deltaResult = Inf;
    
%     %Reset current Time
%     TOF_estimation = pState.initial_tof;
%     pState.tof_current = TOF_estimation;
    figure(2);
    cla;
    %Maximize window
    set(gcf,'WindowState','maximized')

    xlim([initialTime, initialTime + dateSearchSpan])
    ylim([TofLimLow, TofLimHigh])       

    tof_vec = linspace(TofLimLow, TofLimHigh, pSettings.twLabelCount);
    ylims = ylim;
    tof_coords = linspace(ylims(1), ylims(2), pSettings.twLabelCount);

    date_vec = linspace(initialTime, initialTime + dateSearchSpan, pSettings.twLabelCount);
    xlims = xlim;
    date_coords = linspace(xlims(1), xlims(2), pSettings.twLabelCount);

    xticks(date_coords);
    xticklabels(string(secToTime(date_vec, 1, pSettings.baseDate)));

    yticks(tof_coords);
    tickTexts = cell(1, pSettings.twLabelCount);
    for i = 1:pSettings.twLabelCount
        tickTexts(i) =  {sprintf("%.2f y", tof_vec(i) / (365.25 * 86400))};
    end
    yticklabels(tickTexts);

    pState.successfulOrbits = 0;
    pState.testedOrbits = 0;

    pState.twMapInds = [1, 1];

    %-- Solve transfer window using Global Search --
    if pSettings.transferWindowSearchOption == 1
        pSettings.solveDate = 1;
        gs = MultiStart;
        numberOfStartPoints = pSettings.gsPointCount * pSettings.option1pointMultiplier;
        dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
        problem = createOptimProblem('fmincon','x0',[initialTime, pSettings.TOF_average],...
        'objective',dvHandle,'lb',[initialTime, TofLimLow], ...
        'ub',[initialTime + dateSearchSpan, TofLimHigh], 'options', pSettings.opt_dv_global);
        run(gs, problem, numberOfStartPoints);
    end    
    
    %-- Solve transfer window using set dates and optimize TOF --
    if pSettings.transferWindowSearchOption == 2
        for time = linspace(initialTime, initialTime + dateSearchSpan, pSettings.gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong> \n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof_start = linspace(TofLimLow + (TofLimHigh - TofLimLow)/pSettings.option2starts, TofLimHigh, pSettings.option2starts)
                pSettings.solveDate = 0;
                pState.currentTime = time;
                dvHandle = @(tof_in) optimalDVSolver(tof_in, pSettings, []);
                fminsearch(dvHandle, tof_start, pSettings.opt_dv_fminsearch);
            end
        end
    end

    %-- Solve transfer window using a set grid of dates and TOFs --
    if pSettings.transferWindowSearchOption == 3
        pSettings.solveDate = 1;
        for time = linspace(initialTime, initialTime + dateSearchSpan, pSettings.gsPointCount)
            fprintf("Optimizing: <strong>%.1f %%</strong> ---- Best %cV Found: <strong>%.0f m/s</strong>\n", 100 * (time - initialTime)/(dateSearchSpan), 916, deltaResult);
            for tof = linspace(TofLimLow, TofLimHigh, pSettings.gsPointCount)
                pState.currentTime = time;
                pState.tof_current = tof;
                optimalDVSolver([pState.currentTime, pState.tof_current], pSettings, pState.twMapInds);
                pState.twMapInds(2) = pState.twMapInds(2) + 1;
            end
            pState.twMapInds(1) = pState.twMapInds(1) + 1;
            pState.twMapInds(2) = 1;
        end
        %Plot the image
        imagesc([initialTime, initialTime + dateSearchSpan], [TofLimLow, TofLimHigh], pState.twMap);
        %colormap(pSettings.cmap);
    end
    
    fprintf("Optimization Completed! ---- Finalized %cV Found: <strong>%.0f m/s</strong> \n", 916, deltaResult);
    fprintf("Tested %.0f transfer windows out of which %1.f %% (%.0f) were achievable\n", pState.testedOrbits, 100 * (pState.successfulOrbits / pState.testedOrbits), pState.successfulOrbits);
    pState.isRunning = 0;

    toc;
end

%% Hover callback function
function hoverCallback(obj, ~, pSettings, axHandle)
    global pState; %theta_super
    global d_opt paramVector paramVector_opt deltaResult;
    global tfP obP1 obP2_1 obP2_2;

    if pState.isRunning == 0
        pState.isRunning = 1;

        pSettings.plotTransferWindow = 0;

        if pState.axRelative == [0, 0, 0, 0]
            pState.axRelative = get(obj.CurrentAxes, 'Position');
        
            plot(axHandle, pSettings.orbit1(1,:), pSettings.orbit1(2,:), 'LineStyle',':', LineWidth=2);
            plot(axHandle, pSettings.orbit2(1,:), pSettings.orbit2(2,:), 'LineStyle',':', LineWidth=2);

            rectangle(axHandle, 'Position',[-pSettings.bodyRadius, -pSettings.bodyRadius, 2*pSettings.bodyRadius, 2*pSettings.bodyRadius], ...
                     'Curvature',[1 1], 'FaceColor',"yellow")
        end

        %Get the cursor position
        cursorPos = get(obj, 'CurrentPoint');
        windowSize = get(obj, 'Position');
        relativePos = (cursorPos./windowSize(3:4) - pState.axRelative(1:2))./pState.axRelative(3:4);

        xmin = obj.CurrentAxes.XLim(1);
        xSpan = obj.CurrentAxes.XLim(2) - obj.CurrentAxes.XLim(1);
        ymin = obj.CurrentAxes.YLim(1);
        ySpan = obj.CurrentAxes.YLim(2) - obj.CurrentAxes.YLim(1);

        pState.currentTime = xmin + relativePos(1) * xSpan;
        pState.tof_current = ymin + relativePos(2) * ySpan;

        delete(tfP);
        delete(obP1);
        delete(obP2_1);
        delete(obP2_2);

        deltaResult = Inf;
        pSettings.solveDate = 1;
        dv_mouse = optimalDVSolver([pState.currentTime, pState.tof_current], pSettings, []);

        theta_f = paramVector_opt.theta_f;
        r1 = paramVector_opt.r1;
        r2 = paramVector_opt.r2;
        theta1 = paramVector_opt.theta1;
        theta2 = paramVector_opt.theta2;
        nu2_i = paramVector_opt.nu2_i;
        r2_i = paramVector_opt.r2_i;


        if dv_mouse < 1e24
            theta_vec_plot = linspace(pSettings.theta_0, theta_f, pSettings.plotAccuracy);
            theta_super_plot = fGetThetaSuper(theta_vec_plot);
            dT = theta_super_plot(1,2) - theta_super_plot(1,1);

            obP1 = plot(axHandle, cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g');
            obP2_1 = plot(axHandle, cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r');
            obP2_2 = plot(axHandle, cos(pSettings.omega2 + nu2_i) * r2_i, sin(pSettings.omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')   ;   
    
            x = real(cos(theta_vec_plot+theta1) .* fRadiusFunction(d_opt, theta_super_plot, paramVector_opt));
            y = real(sin(theta_vec_plot+theta1) .* fRadiusFunction(d_opt, theta_super_plot, paramVector_opt));
            tfP = plot(axHandle, x, y, "Color", [0.2 0.7 0.2]);


            deltaV_tof = fJerkFunction(d_opt, theta_super_plot, dT, paramVector_opt);

            title(axHandle, sprintf("Transfer date: %s\nUsed TOF: %s\nRequired deltaV: %.0f m/s", secToTime(pState.currentTime, 1, pSettings.baseDate), secToTime(pState.tof_current, 0, []), deltaV_tof));
        
        else
            r1 = paramVector.r1;
            r2 = paramVector.r2;
            theta1 = paramVector.theta1;
            theta2 = paramVector.theta2;
            nu2_i = paramVector.nu2_i;
            r2_i = paramVector.r2_i;

            obP1 = plot(axHandle, cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g');
            obP2_1 = plot(axHandle, cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r');
            obP2_2 = plot(axHandle, cos(pSettings.omega2 + nu2_i) * r2_i, sin(pSettings.omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k');   
            
            title(axHandle, sprintf("Transfer date: %s\nUsed TOF: %s\nNo solution found", secToTime(pState.currentTime, 1, pSettings.baseDate), secToTime(pState.tof_current, 0, [])));
        end

        % Set axes limits
        axHandle.XLim = pSettings.viewportLimits(1:2);
        axHandle.YLim = pSettings.viewportLimits(3:4);
        
        %axis(axHandle, 'square');

        pState.isRunning = 0;
    end
end

function [nu] = nuFromTime(T, n, e, nu_guess)

    M = n*T;
    nu = M;

    if e ~= 0
    
        root = sqrt(-e^3 * (-8 + 24*e - 24*e^2 + 8*e^3 - 9*e*M^2));
    
        nom = -2*e + 2*e^2 + (3*M*e^2 + root)^(2/3);
        denom = e * (3*M*e^2 + root)^(1/3);
    
        E =  nom/denom;
        
        if E <= pi
            nu = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
        else
            nu = 2*pi + 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
        end
    
        if M <= pi
            E2 = 2*atan(tan(M/2)/sqrt((1+e)/(1-e)));
        else
            E2 = 2*pi + 2*atan(tan(M/2)/sqrt((1+e)/(1-e)));
        end


        xnm1 = M;
        fnm1 = E2 - e*sin(E2) - M;
        xn = nu;
        fn = E - e*sin(E) - M;


       
        iter = 0;
        while (abs(fn) > 0.001) && (iter < 100)
            xnp1 = mod(xn - fn * (xn-xnm1) / (fn-fnm1), 2*pi);
            fnm1 = fn;
            xnm1 = xn;
            xn = xnp1;
            fn = nuSolver(xnp1, T, n, e);

            iter = iter + 1;

            if iter == 50
                if abs(fnm1) < abs(fn)
                    fn = nuSolver(nu_guess, T, n, e);
                    xn = nu_guess;
                else
                    fnm1 = nuSolver(nu_guess, T, n, e);
                    xnm1 = nu_guess;            
                end
            end
%             if abs(fn) < abs(minFn)
%                 minFn = fn;
%                 bestXn = xn;
%             end
        end

        nu = xn;%bestXn;

    end

end

%% Solve nu around an orbit at a given time T
function [angleError] = nuSolver(nu_time, T, n, e)
    if nu_time <= pi
        E = 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    else
        E = 2*pi + 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    end

    angleError = E-e*sin(E) - n*T;
end

%% For solving theta at a given time along the transfer orbit
function [timeError] = fThetaFromTime(theta, d, paramVector, targetTime)
    
    
    theta_f = paramVector.theta_f;
    r2 = paramVector.r2;
    
    a = paramVector.a;
    b = paramVector.b;
    c = paramVector.c;
                
    efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
                -tan(paramVector.gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
                paramVector.mju/(r2^4*paramVector.theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
    
    efg = 1/(2*theta_f^6) * paramVector.efg_Mat_1 * efg_Mat_2;
    
    e = efg(1);
    f = efg(2);
    g = efg(3);

    theta_vec = linspace(0, theta, 50);
    dT = theta_vec(2) - theta_vec(1);
    theta_super = [theta_vec; theta_vec.^2; theta_vec.^3; theta_vec.^4; theta_vec.^5; theta_vec.^6];
    
    r = 1 ./ (a + b.*theta_super(1,:) + c.*theta_super(2,:) + d.*theta_super(3,:) + e.*theta_super(4,:) + f.*theta_super(5,:) + g.*theta_super(6,:));

    timeStep_Vec = (r.^2./sqrt(paramVector.mju)) .* sqrt((1./r + 2.*c + 6.*d.*theta_super(1,:) + 12.*e.*theta_super(2,:) + 20.*f.*theta_super(3,:) + 30.*g.*theta_super(4,:))); 
    
    time_t_f = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1));

    timeError = time_t_f - targetTime;

end
