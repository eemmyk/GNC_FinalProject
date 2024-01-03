function [deltaV_o] = optimalDVSolver(inputVec)
    global tof_current currentTime solveDate plot3D theta_vec;
    global d_solution d_minimum d_maximum initial_DeltaV;
    global unfeasibleOrbit;

    if solveDate == 1
        currentTime = inputVec(1);
        tof_current = inputVec(2);
    else
        tof_current = inputVec;
    end

    
    updateParameters(0);

    %Check if the tof is solvable within the limits of the d-coefficient
    t_error_min = trapz(theta_vec, fTimeFunction(0.9*d_maximum, theta_vec, 0)) - tof_current;
    t_error_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0)) - tof_current;
    trueSolution = 1;

    if ((t_error_min < 0) ~= (t_error_max < 0))% && (unfeasibleOrbit == 0)
        opt = optimset('TolFun', 1e2);
        d_solution = fzero(@transferTimeSolution, [d_minimum, d_maximum], opt);

        deltaV_o = trapz(theta_vec, abs(fJerkFunction(d_solution, theta_vec, 0)));
    else
        deltaV_o = 1e24; %A big number
        %fprintf("Unfeasible Time Of Flight Detected\n")
        trueSolution = 0;
    end

    global deltaResult 

    if deltaV_o < deltaResult
        %Get globals from updatedParameters
        global theta2 r2 theta1 r1
        global nu1_i nu2_i r1_i r2_i;
        global theta1_dot theta2_dot;
        global gamma1 gamma2 theta_f;
        
        %Save best results as globals
        global theta2_opt r2_opt;
        global theta1_opt r1_opt;
        global gamma1_opt gamma2_opt
        global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt d_opt;
        global theta1_dot_opt theta2_dot_opt;
        global dateOptimal;


        deltaResult = deltaV_o;
        theta1_opt = theta1;
        theta2_opt = theta2;
        theta_f_opt = theta_f;
        theta1_dot_opt = theta1_dot;
        theta2_dot_opt = theta2_dot;
        gamma1_opt = gamma1;
        gamma2_opt = gamma2;
        d_opt = d_solution;
        r2_opt = r2;
        r1_opt = r1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;

        dateOptimal = currentTime;

%         if solveDate == 1
%             fprintf("Optimal dV: %.0f m/s for transfer date: %.0f s and time: %.0f s\n", deltaV_o, currentTime, tof_current)
%         else
%             fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_current)
%         end
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
% 
%     plotValue = deltaV_o;
%     if ~trueSolution
%         plotValue = 0;
%         color = [0, 0, 0];
%     end
%     

    if plot3D == 1
        plot3(currentTime, tof_current, deltaV_o*trueSolution,'o','Color','k','MarkerSize',10,'MarkerFaceColor', color.*trueSolution)
    else
        plot(tof_current, deltaV_o*trueSolution,'o','Color','k','MarkerSize',10,'MarkerFaceColor', color.*trueSolution)
    end

%     plotValue = deltaV_o;
% 
%     if deltaV_o > 20*initial_DeltaV
%         plotValue = 20*initial_DeltaV;
%     end
% 
%     global contourMap contIndX contIndY contIndLim
%     contourMap(contIndX, contIndY) = log(plotValue);
% 
%     contIndY = contIndY + 1;
%     if contIndY > contIndLim
%         contIndX = contIndX + 1;
%         contIndY = 1;
%     end
end

