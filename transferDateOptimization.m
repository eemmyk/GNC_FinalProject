function [deltaV_o] = transferDateOptimization(inputVec)
    
    global currentTime tof_current d_solution theta_f;
    global d_minimum d_maximum timeResult theta_vec;
 
    currentTime = inputVec(1);
    tof_current = inputVec(2);

    timeResult = Inf;

    updateParameters(0)
    safeTransferAngle = pi;

    if theta_f < safeTransferAngle
        deltaV_o = 1e24; %A big number
        return;
    end

    opt = optimset('TolFun', 1e3);
    d_solution = fzero(@transferTimeOptimization, [d_minimum, d_maximum], opt);

    global thrustFunction thetaDotFunction thetaDotSquareFunction timeFunction radiusFunction;

    syms d theta;

    jerkFunction_nn = @(d_coeff, angle) double(subs(subs(thrustFunction/thetaDotFunction, d, d_coeff), theta, angle));
    
    deltaV_o = trapz(theta_vec, abs(jerkFunction_nn(d_solution, theta_vec)));

    global deltaResult 
    
    if deltaV_o < deltaResult
        %Get globals from updatedParameters
        global theta2 theta2_dot gamma2 r2 P2 theta1 theta1_dot gamma1 r1 P1
        global nu1_i nu2_i r1_i r2_i;
    
        %Save best results as globals
        global theta2_opt r2_opt tof_optimal theta2_dot_opt gamma2_opt P2_opt;
        global theta1_opt r1_opt theta1_dot_opt gamma1_opt P1_opt theta_f_opt;
        global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt dateOptimal d_opt;
    
        global thrustFunction_opt thetaDotFunction_opt thetaDotSquareFunction_opt;
        global radiusFunction_opt timeFunction_opt;

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

        fprintf("Optimal dV: %.0f m/s -- transfer date: %.0f s -- time of flight: %.0f s\n", deltaV_o, inputVec(1), inputVec(2))
    end

    plot3(currentTime, tof_current, deltaV_o,'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF')
    pause(0.001)
  
end

