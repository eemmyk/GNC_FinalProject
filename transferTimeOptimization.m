function [t_error] = transferTimeOptimization(d_in)
    global theta_f theta_0 tf timeFunction  integralApproximationSteps;
    syms d theta;

    timeFunction_n = subs(timeFunction, d, d_in);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));

    %Transfer Time
    theta_vec = linspace(theta_0, theta_f, integralApproximationSteps);
    timeCurve = zeros(2, integralApproximationSteps);
    for i = 1:integralApproximationSteps
        timeCurve(:,i) = [theta_vec(i); timeFunction_nn(theta_vec(i))];
    end
    
    %Getting rid of the imaginary numbers let's the code run, but the
    %solutions are not possible with only retrograde/prograde thrust
    %time_t = abs(trapz(timeCurve(1, :), timeCurve(2,:)));
    
    time_t = real(trapz(timeCurve(1, :), timeCurve(2,:)));

    t_error = time_t - tf;

    fprintf("Guessed d: %e, Remaining TOF error: %.0f\n", d_in, t_error);

end

