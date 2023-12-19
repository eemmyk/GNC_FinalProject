function [deltaV_o] = deltaVOptimization(d_in)
    global theta_f theta_0 thrustFunction thetaDotFunction integralApproximationSteps;
    
    syms d theta;

    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_in);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    thrustFunction_n = subs(thrustFunction, d, d_in);
    thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));


    theta_vec = linspace(theta_0, theta_f, integralApproximationSteps);
    jerk = zeros(2, integralApproximationSteps);
    thrust = zeros(2, integralApproximationSteps);
    

    for i = 1:integralApproximationSteps
        jerk(:,i) = [theta_vec(i); jerkFunction_nn(theta_vec(i))];
        thrust(:,i) = [theta_vec(i); thrustFunction_nn(theta_vec(i))];
    end

    deltaV_o = trapz(jerk(1, :), abs(jerk(2,:)));
    
    fprintf("Guessed d: %e, Solution dV: %.0f\n", d_in, deltaV_o);

end

