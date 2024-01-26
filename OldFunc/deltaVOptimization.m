function [deltaV_o] = deltaVOptimization(d_in)
    global theta_f theta_0 thrustFunction thetaDotFunction intApprox;
    
    syms d theta;

    jerkFunction_n = subs(thrustFunction/thetaDotFunction, d, d_in);
    jerkFunction_nn = @(angle) double(subs(jerkFunction_n, theta, angle));
    %thrustFunction_n = subs(thrustFunction, d, d_in);
    %thrustFunction_nn = @(angle) double(subs(thrustFunction_n, theta, angle));


    theta_vec = linspace(theta_0, theta_f, intApprox);
%     jerk = zeros(2, intApprox);
%     thrust = zeros(2, intApprox);
%     
% 
%     for i = 1:intApprox
%         jerk(:,i) = [theta_vec(i); jerkFunction_nn(theta_vec(i))];
%         thrust(:,i) = [theta_vec(i); thrustFunction_nn(theta_vec(i))];
%     end
    jerk = [theta_vec; jerkFunction_nn(theta_vec)];
    %thrust = [theta_vec; thrustFunction_nn(theta_vec)];

    
    deltaV_o = trapz(jerk(1, :), abs(jerk(2,:)));
    
    %Save best result into a global variable
    global interDeltaResult;
    if deltaV_o < interDeltaResult
        interDeltaResult = deltaV_o;
        %fprintf("Guessed d: %e, Solution dV: %.0f\n", d_in, deltaV_o);
    end
    
    fprintf("Guessed d: %e, Solution dV: %.0f\n", d_in, deltaV_o);
    
    
end

