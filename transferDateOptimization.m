function [cost_o] = transferDateOptimization(inputVec)
    
    global currentTime tof_current d_solution theta_0 theta_f intApprox costResult timeResult;

    currentTime = inputVec(1);
    tof_current = inputVec(2);
    %d_solution = inputVec(3);

    timeResult = Inf;

    updateParameters(0)

    opt = optimset('TolFun', 1e-9, 'TolX', 1e-9);
    d_solution = fzero(@transferTimeOptimization, d_solution, opt);

    
%     global timeFunction;
%     syms d theta;
% 
%     timeFunction_n = subs(timeFunction, d, d_solution);
%     timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));
% 
%     theta_vec = linspace(theta_0, theta_f, intApprox);
%     timeCurve = [theta_vec; timeFunction_nn(theta_vec)];
%     time_t = real(trapz(timeCurve(1, :), timeCurve(2,:)));
%     t_error = time_t - tof_current;

    global thrustFunction thetaDotFunction thetaDotSquareFunction timeFunction radiusFunction;

    syms d theta;

    theta_vec = linspace(theta_0, theta_f, intApprox);
    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_solution);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    jerk = [theta_vec; jerkFunction_nn(theta_vec)];
    
    deltaV_o = trapz(jerk(1, :), abs(jerk(2,:)));

    %Cost functions
    %w1 = 50; %TOF error cost
    %w2 = 0.000026; %DeltaV cost
    %cost_o = w1*timeResult^2 + w2*deltaV_o^2;

    cost_o = deltaV_o;

    %Get globals from updatedParameters
    global theta2 theta2_dot gamma2 r2 P2 theta1 theta1_dot gamma1 r1 P1
    global nu1_i nu2_i r1_i r2_i;

    %Save best results as globals
    global deltaResult theta2_opt r2_opt tof_optimal theta2_dot_opt gamma2_opt P2_opt;
    global theta1_opt r1_opt theta1_dot_opt gamma1_opt P1_opt theta_f_opt;
    global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt dateOptimal d_opt;

    global thrustFunction_opt thetaDotFunction_opt thetaDotSquareFunction_opt;
    global radiusFunction_opt timeFunction_opt;

    if cost_o < costResult
        costResult = cost_o;
        deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta_f_opt = theta_f;
        theta2_dot_opt = theta2_dot;
        theta1_dot_opt = theta1_dot;
        gamma2_opt = gamma2;
        gamma1_opt = gamma1;
        d_opt = d_solution;
        tof_optimal = tof_current;
        dateOptimal = currentTime;
        r2_opt = r2;
        r1_opt = r1;
        P2_opt = P2;
        P1_opt = P1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;

        thrustFunction_opt =  thrustFunction;
        thetaDotFunction_opt = thetaDotFunction; 
        thetaDotSquareFunction_opt = thetaDotSquareFunction; 
        timeFunction_opt = timeFunction;
        radiusFunction_opt = radiusFunction;

        %fprintf("Optimal dV: %.0f m/s -- transfer date: %.0f s -- time of flight: %.0f s\n", deltaV_o, inputVec(1), inputVec(2))
    end

    plot3(currentTime, tof_current, deltaV_o);
    pause(0.001)
    
    fprintf("Optimal dV: %.0f m/s -- transfer date: %.0f s -- time of flight: %.0f s\n", deltaV_o, inputVec(1), inputVec(2))
   
    %fprintf("Optimal dV: %.0f m/s for transfer date: %.0f s\n", deltaV_o, currentTime)
    %printf("Total cost: %e -- Date shifted: %e -- Time of Flight: %e\n", cost_o, currentTime, tof_current)
    
end

