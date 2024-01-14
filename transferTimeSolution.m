function [t_error] = transferTimeSolution(d_in, paramVector, tof_current, theta_super) 
    %Own trapz implementation
    dT = theta_super(1,2) - theta_super(1,1);
    timeStep_Vec = fTimeFunction(d_in, theta_super, paramVector);
    time_t = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1));
    
    %Some solutions are still imaginary for some values of d.
    %So far values always converge into real numbers after a while.
    if isreal(time_t)
        t_error = time_t - tof_current;
    else
        t_error = sign(real(time_t - tof_current))*1e24; %Arbitary large value with same sign
    end
end

