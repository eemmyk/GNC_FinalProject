function [deltaV_o] = optimalDVSolver(inputVec, pSettings)
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

    for N_current = N_start:N_end
        
        pState.N = N_current;
        [resultVector, paramVector] = updateParameters(0, pSettings);
        
        d_minimum = resultVector(1);
        d_maximum = resultVector(2);

        try
            tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, pState.tof_current, theta_super);
            d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], pSettings.opt_tof_fzero);
        
            dT = theta_super(1,2) - theta_super(1,1);
            deltaV_o_Vec = abs(fJerkFunction(d_solution, theta_super, paramVector));
            deltaV_o = dT * (deltaV_o_Vec(1) + deltaV_o_Vec(end)) / 2 + dT * sum(deltaV_o_Vec(2:end-1));


            %deltaV_o = trapz(theta_super(1,:), abs(fJerkFunction(d_solution, theta_super, paramVector)));
            trueSolution = 1;
        catch
            d_solution = 0;
            deltaV_o = 1e24; %A big number
        end

        if deltaV_o < localBestDV
            localBestDV = deltaV_o;
        end

        if deltaV_o < deltaResult
            %Save best results as globals
            global dateOptimal tof_optimal d_opt;
            global paramVector_opt
    
            deltaResult = deltaV_o;
            d_opt = d_solution;
    
            dateOptimal = pState.currentTime;
            tof_optimal = pState.tof_current;
    
            paramVector_opt = paramVector;
        end
    end

    if pSettings.plotTransferWindow == 1
        colorScaleDown = 4;

        R_Multiplier = (3*(localBestDV/pState.initial_DeltaV)^2 - 2*(localBestDV/pState.initial_DeltaV)^3);
        G_Multiplier = 1-(3*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^2 - 2*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^3);
    
        if localBestDV > (1+colorScaleDown)*pState.initial_DeltaV
            color = [1 0 0];
        elseif localBestDV > pState.initial_DeltaV
            color = [1 G_Multiplier^4 0];
        elseif localBestDV <= pState.initial_DeltaV
            color = [R_Multiplier^4 1 0];
        end

        color = color.*(0.75 + 0.25 * trueSolution);

        rectangle('Position',[pState.currentTime-0.5*pSettings.tfWindowPixelsX, pState.tof_current-0.5*pSettings.tfWindowPixelsY, ...
                              pSettings.tfWindowPixelsX, pSettings.tfWindowPixelsY], 'FaceColor', color, 'EdgeColor', color);
                           
    end   

    %Reset value
    pState.N = N_start;
    
    if trueSolution == 0
        pState.failedOrbits = pState.failedOrbits + 1;
    end
end