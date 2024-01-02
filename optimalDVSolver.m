function [deltaV_o] = optimalDVSolver(tof_in)
    global tof_current  theta_f currentTime plotDV_3D theta_vec;
    global d_solution d_minimum d_maximum initial_DeltaV timeResult;
    
    tof_current = tof_in;
    timeResult = Inf;
    
%     global N
%     N = 0;

    updateParameters(0);

%     safeTransferAngle = pi;

%     if theta_f < safeTransferAngle
%         deltaV_o = 1e24; %A big number
%         return;
%     end


%     timeFunction
%     d_minimum
%     d_maximum
%     theta1
%     theta2
%     theta_f
%     tof_current

    %Check if the tof is solvable within the limits of the d-coefficient
    t_error_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec)) - tof_current;
    t_error_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec)) - tof_current;
    trueSolution = 1;

    if (t_error_min < 0) ~= (t_error_max < 0)
        opt = optimset('TolFun', 1e2);
        d_solution = fzero(@transferTimeOptimization, [d_minimum, d_maximum], opt);

        deltaV_o = trapz(theta_vec, abs(fJerkFunction(d_solution, theta_vec)));
    else
        deltaV_o = 1e24; %A big number
        %fprintf("Unfeasible Time Of Flight Detected\n")
        trueSolution = 0;
    end

    global deltaResult 

    if deltaV_o < deltaResult
        %Get globals from updatedParameters
        global radiusFunction;
        global theta2 r2 theta1 r1
        global nu1_i nu2_i r1_i r2_i;
        
        %Save best results as globals
        global theta2_opt r2_opt tof_optimal;
        global theta1_opt r1_opt;
        global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt d_opt;
        global radiusFunction_opt dateOptimal;


        deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta_f_opt = theta_f;
        d_opt = d_solution;
        tof_optimal = tof_in;
        r2_opt = r2;
        r1_opt = r1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;


        radiusFunction_opt = radiusFunction;

        dateOptimal = currentTime;

        fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)
    end

    R_Multiplier = (3*(deltaV_o/initial_DeltaV)^2 - 2*(deltaV_o/initial_DeltaV)^3);
    G_Multiplier = 1-(3*((deltaV_o-initial_DeltaV)/initial_DeltaV)^2 - 2*((deltaV_o-initial_DeltaV)/initial_DeltaV)^3);
    
    if deltaV_o > 2*initial_DeltaV
        color = [1 0 0];
    elseif deltaV_o > initial_DeltaV
        color = [1 G_Multiplier 0];
    elseif deltaV_o <= initial_DeltaV
        color = [R_Multiplier 1 0];
    end

    plotValue = deltaV_o;
    if ~trueSolution
        plotValue = 0;
        color = [0, 0, 0];
    end
    
    if plotDV_3D == 1
        plot3(currentTime, tof_current, plotValue,'o','Color','k','MarkerSize',10,'MarkerFaceColor', color)
    else
        plot(tof_current, plotValue,'o','Color','k','MarkerSize',10,'MarkerFaceColor', color)
    end
    pause(0.001);
end

