function [deltaV_o] = optimalDVSolver(tof_in)
    
    %fprintf("transfer time: %.0f s\n", tof_in)
    global tof_current currentTime timeResult;

    tof_current = tof_in;
    
    updateParameters()
    
    %Initial guess for d coefficient:
    global d_solution theta_0 theta_f intApprox;
    %d_guess = d_solution;

    % Use function to find best dV value for given tf

    %d_fuelOptimal_guess = d_guess;
    
    %timeResult = Inf;

    opt = optimset('TolFun',1e1);
    %global interDeltaResult

    %d_fuelOptimal = fminsearch(@deltaVOptimization, d_fuelOptimal_guess, opt);
    d_solution = fzero(@transferTimeOptimization, d_solution, opt);

    %jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_fuelOptimal);
    %jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    %thrustFunction_n = subs(thrustFunction, d, d_fuelOptimal);
    %thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));

    %theta_vec = linspace(theta_0, theta_f, intApprox);
%     jerk = zeros(2, intApprox);
%     thrust = zeros(2, intApprox);
%     for i = 1:intApprox
%         jerk(:,i) = [theta_vec(i); double(subs(jerkFunction_n, theta, theta_vec(i)))];
%         thrust(:,i) = [theta_vec(i); double(subs(thrustFunction_n, theta, theta_vec(i)))];
%     end

    %jerk = [theta_vec; jerkFunction_nn(theta_vec)];
    %thrust = [theta_vec; thrustFunction_nn(theta_vec)];    

    %deltaV_o = interDeltaResult; %trapz(jerk(1, :), abs(jerk(2,:)));

    global thrustFunction thetaDotFunction;

    syms d theta;

    theta_vec = linspace(theta_0, theta_f, intApprox);
    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_solution);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    jerk = [theta_vec; jerkFunction_nn(theta_vec)];
    %thrust = [theta_vec; thrustFunction_nn(theta_vec)];

    
    deltaV_o = trapz(jerk(1, :), abs(jerk(2,:)));

    %Get globals from updatedParameters
    global theta2 theta2_dot gamma2 r2 P2 theta1 theta1_dot gamma1 r1 P1
    global nu1_i nu2_i r1_i r2_i;

    %Save best results as globals
    global deltaResult theta2_opt r2_opt tof_optimal theta2_dot_opt gamma2_opt P2_opt;
    global theta1_opt r1_opt theta1_dot_opt gamma1_opt P1_opt;
    global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt d_fuelOptimal;

    if deltaV_o < deltaResult
        deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta2_dot_opt = theta2_dot;
        theta1_dot_opt = theta1_dot;
        gamma2_opt = gamma2;
        gamma1_opt = gamma1;
        d_fuelOptimal = d_solution;
        tof_optimal = tof_in;
        r2_opt = r2;
        r1_opt = r1;
        P2_opt = P2;
        P1_opt = P1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;


        %fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)
    end
    fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)

end

