function [deltaV_o] = optimalDVSolver(tof_in)
    
    %fprintf("transfer time: %.0f s\n", tof_in)
    global tof_current N theta_f currentTime;

    tof_current = tof_in;
    
    %N = 0;

    updateParameters(0)
    safeTransferAngle = pi;

    if theta_f < safeTransferAngle
        %N = N+1;
        %updateParameters(1);
        deltaV_o = 1e24;
        return;
    end
    
    %tof_current = tof_in;


    %Initial guess for d coefficient:
    global d_solution theta_0 theta_f intApprox d_minimum d_maximum;
    %d_guess = d_solution;

    % Use function to find best dV value for given tf

    %d_fuelOptimal_guess = d_guess;
    
    %timeResult = Inf;

    opt = optimset('TolFun', 1e3);
    %global interDeltaResult

    %d_fuelOptimal = fminsearch(@deltaVOptimization, d_fuelOptimal_guess, opt);
    d_solution = fzero(@transferTimeOptimization, [d_minimum, d_maximum], opt);
    %d_out = fminsearch(@transferTimeOptimization, d_solution);

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

    global thrustFunction thetaDotFunction thetaDotSquareFunction;
    global timeFunction radiusFunction;

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
    global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt d_opt;
    global thrustFunction_opt thetaDotFunction_opt thetaDotSquareFunction_opt;
    global radiusFunction_opt timeFunction_opt;

    if deltaV_o < deltaResult
        deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta_f_opt = theta_f;
        theta2_dot_opt = theta2_dot;
        theta1_dot_opt = theta1_dot;
        gamma2_opt = gamma2;
        gamma1_opt = gamma1;
        d_opt = d_solution;
        tof_optimal = tof_in;
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


        fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)
    end
    %fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)
    
    plot3(currentTime, tof_current, deltaV_o,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    %plot(tof_in, deltaV_o,'or', 'MarkerSize',2,'MarkerFaceColor','r')
    pause(0.001);
end

