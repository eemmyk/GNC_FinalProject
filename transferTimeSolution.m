function [t_error] = transferTimeSolution(d_in, paramVector, tof_current, theta_vec)    
    time_t = trapz(theta_vec, fTimeFunction(d_in, theta_vec, paramVector));
    

    %Some solutions are still imaginary for some values of d.
    %So far values always converge into real numbers after a while.
    if isreal(time_t)
        t_error = time_t - tof_current;
    else
        t_error = sign(real(time_t - tof_current))*1e24; %Arbitary large value with same sign
    end
end

