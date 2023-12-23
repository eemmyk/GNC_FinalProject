function [t_error] = transferTimeOptimization(d_in)
    global theta_f theta_0 tf timeFunction  intApprox N tof_current d_solution;
    
    syms d theta;

    timeFunction_n = subs(timeFunction, d, d_in);
    timeFunction_nn = @(angle) double(subs(timeFunction_n, theta, angle));

    %Transfer Time
    theta_vec = linspace(theta_0, theta_f, intApprox);
%     timeCurve = zeros(2, intApprox);
%     for i = 1:intApprox
%         timeCurve(:,i) = [theta_vec(i); timeFunction_nn(theta_vec(i))];
%     end
%     
    timeCurve = [theta_vec; timeFunction_nn(theta_vec)];

    %Getting rid of the imaginary numbers let's the code run, but the
    %solutions are not possible with only retrograde/prograde thrust
    %time_t = abs(trapz(timeCurve(1, :), timeCurve(2,:)));
    
%     if ~isreal(timeCurve)
%         N = mod(N + 1, 3);
%         d_solution = 1e-9;
%         updateParameters()
        time_t = real(trapz(timeCurve(1, :), timeCurve(2,:)));
        t_error = time_t - tof_current;
%     else
%         %time_t = real(trapz(timeCurve(1, :), timeCurve(2,:)));
%         time_t = trapz(timeCurve(1, :), timeCurve(2,:));
%         t_error = time_t - tof_current;
%     end

    %Save best result into a global variable
    global timeResult;
    if abs(t_error) < abs(timeResult)
        timeResult = t_error;
        %fprintf("Guessed d: %e, Remaining TOF error: %.0f\n", d_in, t_error);
    end

    fprintf("Guessed d: %e, Remaining TOF error: %.0f\n", d_in, t_error);
    
    
end

