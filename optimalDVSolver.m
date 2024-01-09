function [deltaV_o] = optimalDVSolver(inputVec)
    global tof_current currentTime solveDate plotTransferWindow theta_vec;
    global d_solution d_minimum d_maximum initial_DeltaV;
    global opt_tof_fzero;
    global deltaResult 
    global tfWindowPixelsX tfWindowPixelsY;

    %global tw_graph_ind tw_graph unfeasibleOrbit;

    global paramVector

    if solveDate == 1
        currentTime = inputVec(1);
        tof_current = inputVec(2);
    else
        tof_current = inputVec;
    end

    updateParameters(0);


%     %Check if the tof is solvable within the limits of the d-coefficient
%     t_error_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, paramVector)) - tof_current;
%     t_error_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, paramVector)) - tof_current;

    %trueSolution = 0;
    
    %Block Impossible orbits that escape to infinite, pass throught the
    %planet or require lateral thrust (maybe on the last one)
%     realSolutions = isreal(t_error_min) && isreal(t_error_max);
%     zeroCrossingSolutions = ((t_error_min < 0) ~= (t_error_max < 0));

%     if realSolutions && zeroCrossingSolutions % && (unfeasibleOrbit == 0)
%         tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, tof_current, theta_vec);
%         d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], opt_tof_fzero);
% 
%         deltaV_o = trapz(theta_vec, abs(fJerkFunction(d_solution, theta_vec, paramVector)));
%     else
%         deltaV_o = 1e24; %A big number
%         %fprintf("Unfeasible Time Of Flight Detected\n")
%         trueSolution = 0;
%     end
%     
    %while (trueSolution == 0) && (N < 4)
        try
            tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, tof_current, theta_vec);
            d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], opt_tof_fzero);
    
            deltaV_o = trapz(theta_vec, abs(fJerkFunction(d_solution, theta_vec, paramVector)));
            trueSolution = 1;
        catch
            %N = N+1;
            %updateParameters(0);
%             t_min = fTimeFunction(d_maximum, theta_vec, paramVector)
%             t_max = fTimeFunction(d_minimum, theta_vec, paramVector)
%             t_error_min = trapz(theta_vec, t_min) - tof_current
%             t_error_max = trapz(theta_vec, t_max) - tof_current
%             ang = [max(t_min - real(t_min)), min(t_min - real(t_min))]
%             ang2 = [max(t_max - real(t_max)), min(t_max - real(t_max))]
            deltaV_o = 1e24; %A big number
            trueSolution = 0;
        end
    %end
    %N = 0;
 
    if deltaV_o < deltaResult
        %Save best results as globals
        global dateOptimal tof_optimal d_opt;
        global paramVector_opt

        deltaResult = deltaV_o;
        d_opt = d_solution;

        dateOptimal = currentTime;
        tof_optimal = tof_current;

        paramVector_opt = paramVector;
        
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

    if plotTransferWindow == 1
        %plot3(currentTime, tof_current, deltaV_o*trueSolution,'o','Color','k','MarkerSize',6,'MarkerFaceColor', color.*trueSolution)
        rectangle('Position',[currentTime-0.5*tfWindowPixelsX, tof_current-0.5*tfWindowPixelsY, tfWindowPixelsX, tfWindowPixelsY], ...
                  'FaceColor', color.*trueSolution, 'EdgeColor',color.*trueSolution)
        %pause(0);
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

