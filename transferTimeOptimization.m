function [t_error] = transferTimeOptimization(d_in)
    global tof_current theta_vec;
    global timeFunction_nn

    timeCurve = [theta_vec; timeFunction_nn(d_in, theta_vec)];

    time_t = trapz(timeCurve(1, :), timeCurve(2,:));
    t_error = time_t - tof_current;

    %Save best result into a global variable
    global timeResult;
    if abs(t_error) < abs(timeResult)
        timeResult = t_error;
        fprintf("Guessed d: %e, Remaining TOF error: %.0f\n", d_in, t_error);
    end    
end

