function [t_error] = trueTofSolver(currentTime_in)
    global currentTime;
    
    currentTime = currentTime_in;
    updateParameters()
    
    %Initial guess for d coefficient:
    global d_solution theta_0 theta_f intApprox timeFunction tof_current;
    
    opt = optimset('TolFun',1e1);
    
    d_solution = fzero(@transferTimeOptimization, d_solution, opt);

    syms d theta;
    timeFunction_n = subs(timeFunction, d, d_solution);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));

    %Transfer Time
    theta_vec = linspace(theta_0, theta_f, intApprox);
    timeCurve = [theta_vec; timeFunction_nn(theta_vec)];
    time_t = trapz(timeCurve(1, :), timeCurve(2,:));
    t_error = time_t - tof_current;
    
    %Get globals from updatedParameters
    global theta2 theta2_dot gamma2 r2 P2 theta1 theta1_dot gamma1 r1 P1
    global nu1_i nu2_i r1_i r2_i;

    %Save best results as globals
    global deltaResult theta2_opt r2_opt timeOptimal theta2_dot_opt gamma2_opt P2_opt;
    global theta1_opt r1_opt theta1_dot_opt gamma1_opt P1_opt;
    global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt;

    if deltaV_o < deltaResult
        %deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta2_dot_opt = theta2_dot;
        theta1_dot_opt = theta1_dot;
        gamma2_opt = gamma2;
        gamma1_opt = gamma1;
        %d_solution = d_fuelOptimal;
        timeOptimal = currentTime_in;
        r2_opt = r2;
        r1_opt = r1;
        P2_opt = P2;
        P1_opt = P1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;


        fprintf("Transfer window solution found at %e ", timeOptimal)
    end
    %fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)

end

