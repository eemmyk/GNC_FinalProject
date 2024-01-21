function [deltaV_o] = optimalDVSolver(inputVec, pSettings, twMapInds)
    global d_solution deltaResult theta_super pState; 

    if pSettings.solveDate == 1
        pState.currentTime = inputVec(1);
        pState.tof_current = inputVec(2);
    else
        pState.tof_current = inputVec;
    end

    trueSolution = 0;
    N_start = pState.N;
    localBestDV = Inf;

    %Set limits for multiorbit search
    if pSettings.useMultiorbitFilling == 1
        N_end = pSettings.maxDepthN;
    else
        N_end = N_start;
    end

    pState.testedOrbits = pState.testedOrbits + 1;
    
    d_solution = 0;
    deltaV_o = 1e24; %A big number

    for N_current = N_start:N_end
        
        pState.N = N_current;
        [resultVector, paramVector] = updateParameters(0, pSettings);
        dT = theta_super(1,2) - theta_super(1,1);
        tof_current = pState.tof_current;

        d_minimum = resultVector(1);
        d_maximum = resultVector(2);
        realOrbit = resultVector(4);

         if ~realOrbit
            continue
         end
    
        x_min = d_minimum;
        timeStep_Vec = fTimeFunction(x_min, theta_super, paramVector);
        f_min = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
        x_max = d_maximum;
        timeStep_Vec = fTimeFunction(x_max, theta_super, paramVector);
        f_max = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
        

        %%%
    
%         figure(4)
%         hold on
% 
%         r1 = paramVector.r1;
%         r2 = paramVector.r2;
%         theta_f = paramVector.theta_f;
%         theta1 = paramVector.theta1;
%         theta2 = paramVector.theta2;
%         nu2_i = paramVector.nu2_i;
%         r2_i = paramVector.r2_i;
% 
%         plot(pSettings.orbit1(1,:), pSettings.orbit1(2,:), 'LineStyle',':', LineWidth=2);
%         plot(pSettings.orbit2(1,:), pSettings.orbit2(2,:), 'LineStyle',':', LineWidth=2);
%         
%         plot(cos(theta1) * r1, sin(theta1) * r1,'or', 'MarkerSize',5,'MarkerFaceColor','g')
%         plot(cos(theta2) * r2, sin(theta2) * r2,'or', 'MarkerSize',5,'MarkerFaceColor','r')
%         plot(cos(pSettings.omega2 + nu2_i) * r2_i, sin(pSettings.omega2 + nu2_i) * r2_i,'or', 'MarkerSize',5,'MarkerFaceColor','k')
% 
%         theta_vec_plot = linspace(0, theta_f, 100);
%         %theta_vec_super = fGetThetaSuper(theta_vec_plot);
%     
% 
%         fprintf("Inital Time of Flight guess not achievable\n")
% 
%         x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_minimum*1.01, theta_super, paramVector);
%         y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_minimum*1.01, theta_super, paramVector);
%         time_max = trapz(theta_vec_plot, fTimeFunction(d_minimum*1.01, theta_super, paramVector)) - tof_current;
%         plot(x, y, "Color", [0.5 0.9 0.5]);
% 
%         x = cos(theta_vec_plot+theta1) .* fRadiusFunction(d_maximum, theta_super, paramVector);
%         y = sin(theta_vec_plot+theta1) .* fRadiusFunction(d_maximum, theta_super, paramVector);
%         time_min = trapz(theta_vec_plot, fTimeFunction(d_maximum, theta_super, paramVector)) - tof_current;
%         plot(x, y, "Color", [0.5 0.9 0.5]);

%%%%%


        crossing = (f_min < 0) ~= (f_max < 0);
        if ~crossing % || imaginary
            continue
        end
        
        d_minimum_in = d_minimum;
        d_maximum_in = d_maximum;
        
%         if (d_minimum < 0) ~= (d_maximum < 0)
%             timeStep_Vec = fTimeFunction(0, theta_super, paramVector);
%             f_zero = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
%             if f_zero > 0
%                 d_minimum_in = 0;
%                 f_min = f_zero;
%             else
%                 d_maximum_in = 0;
%                 f_max = f_zero;
%             end
%         end

%         crossing = (f_min < 0) ~= (f_max < 0);
%         if ~crossing % || imaginary
%             continue
%         end


        tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, pState.tof_current, theta_super);
        d_solution = fMyFastZero(tfTimeHandle, [d_minimum_in, d_maximum_in], pSettings.opt_tof_fzero);

        dT = theta_super(1,2) - theta_super(1,1);
        deltaV_o_Vec = abs(fJerkFunction(d_solution, theta_super, paramVector));
        deltaV_o = dT * (deltaV_o_Vec(1) + deltaV_o_Vec(end)) / 2 + dT * sum(deltaV_o_Vec(2:end-1));
        trueSolution = 1;

        if deltaV_o <= localBestDV
            localBestDV = deltaV_o;
        else
            break
        end

        if deltaV_o < deltaResult
            %Save best results as globals
            global dateOptimal tof_optimal d_opt;
            global paramVector_opt resultVector_opt
    
            deltaResult = deltaV_o;
            d_opt = d_solution;
    
            dateOptimal = pState.currentTime;
            tof_optimal = pState.tof_current;
    
            paramVector_opt = paramVector;
            resultVector_opt = resultVector;
        end

    end

    if pSettings.plotTransferWindow == 1
        
        if pSettings.transferWindowSearchOption ~= 3
            colorScaleDown = 4;
    
            R_Multiplier = (3*(localBestDV/pState.initial_DeltaV)^2 - 2*(localBestDV/pState.initial_DeltaV)^3);
            G_Multiplier = 1-(3*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^2 - 2*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^3);
        
            if localBestDV > (1+colorScaleDown)*pState.initial_DeltaV
                color = [1 0 0];
            elseif localBestDV > pState.initial_DeltaV
                color = [1 sqrt(G_Multiplier) 0];
            elseif localBestDV <= pState.initial_DeltaV
                color = [R_Multiplier^2 1 0];
            end
    
            color = color.*(0.75 + 0.25 * trueSolution);
            rectangle('Position',[pState.currentTime-0.5*pSettings.tfWindowPixelsX, pState.tof_current-0.5*pSettings.tfWindowPixelsY, ...
                       pSettings.tfWindowPixelsX, pSettings.tfWindowPixelsY], 'FaceColor', color, 'EdgeColor', color);
        else
            if localBestDV > 2*pState.initial_DeltaV
                value = 0.01;
            else
                value = 0.01 + 1 - (localBestDV/(2*pState.initial_DeltaV));
            end
            
            pState.twMap(twMapInds(2), twMapInds(1)) = value * trueSolution;

        end
    
    end   

    %Reset value
    pState.N = N_start;
    
    if trueSolution == 0
        pState.failedOrbits = pState.failedOrbits + 1;
    end
end

% function root = myFzero(func, x0, x1, tol, max_iter)
%     global theta_super pState
% 
%     dT = theta_super(1,2) - theta_super(1,1);
% 
% 
%     if nargin < 4
%         tol = 1e-6;
%     end
%     
%     if nargin < 5
%         max_iter = 100;
%     end
%     
%     for iteration = 1:max_iter
%         f0_vec = func(x0);
%         f0 = dT * (f0_vec(1) + f0_vec(end)) / 2 + dT * sum(f0_vec(2:end-1)) - pState.tof_current;
% 
%         f1_vec = func(x1);
%         f1 = dT * (f1_vec(1) + f1_vec(end)) / 2 + dT * sum(f1_vec(2:end-1)) - pState.tof_current;
%         
%         if abs(f1) < tol
%             root = x1;
%             return;
%         end
%         
%         % Check for sign change and apply bisection
%         if (f0 < 0) ~= (f1 < 0)
%             x2 = (x0 + x1) / 2;
%         else
%             % Apply secant method
%             x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
%             
%             % Check for convergence
% %             if abs(x2 - x1) < tol
% %                 root = x2;
% %                 return;
% %             end
%             
%             % Apply inverse quadratic interpolation
%             x3 = x0 - (f0 * (x1 - x0)^2) / (f1 - 2*f0 + f1);
% %             if abs(x3 - x2) < tol
% %                 root = x3;
% %                 return;
% %             end
% 
%             f3_vec = func(x3);
%             f3 = dT * (f3_vec(1) + f3_vec(end)) / 2 + dT * sum(f3_vec(2:end-1)) - pState.tof_current;
%             
%             % Choose the best approximation based on the bisection method
%             if (f1 < 0) ~= (f3 < 0)
%                 x2 = (x1 + x3) / 2;
%             end
%         end
%         
%         % Update iterates
%         x0 = x1;
%         x1 = x2;
%     end
%     
%     error('Root not found within the maximum number of iterations.');
% end