function [deltaV_o] = optimalDVSolver(inputVec, pSettings)
    global d_solution paramVector resultVector deltaResult theta_vec pState; 

    global faultCount

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
        updateParameters(0, pSettings);
        
        d_minimum = resultVector(1);
        d_maximum = resultVector(2);
    
        try
            tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, pState.tof_current, theta_vec);
            d_solution = fzero(tfTimeHandle, [d_minimum, d_maximum], pSettings.opt_tof_fzero);
    
            deltaV_o = trapz(theta_vec, abs(fJerkFunction(d_solution, theta_vec, paramVector)));
            trueSolution = 1;
        catch
            d_solution = 0;
            deltaV_o = 1e24; %A big number
        end

        if deltaV_o < localBestDV
            localBestDV = deltaV_o;

            if pSettings.plotTransferWindow == 1
                R_Multiplier = (3*(deltaV_o/pState.initial_DeltaV)^2 - 2*(deltaV_o/pState.initial_DeltaV)^3);
                G_Multiplier = 1-(3*((deltaV_o-pState.initial_DeltaV)/pState.initial_DeltaV)^2 - 2*((deltaV_o-pState.initial_DeltaV)/pState.initial_DeltaV)^3);
            
                if deltaV_o > 2*pState.initial_DeltaV
                    color = [1 0 0];
                elseif deltaV_o > pState.initial_DeltaV
                    color = [1 G_Multiplier 0];
                elseif deltaV_o <= pState.initial_DeltaV
                    color = [R_Multiplier 1 0];
                end
        
                rectangle('Position',[pState.currentTime-0.5*pSettings.tfWindowPixelsX, pState.tof_current-0.5*pSettings.tfWindowPixelsY, ...
                                      pSettings.tfWindowPixelsX, pSettings.tfWindowPixelsY], 'FaceColor', color.*(0.5 + 0.5 * trueSolution), 'EdgeColor',color.*(0.5 + 0.5 * trueSolution));
                                   
            end   
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

    %Reset value
    pState.N = N_start;
    
    if trueSolution == 0
        pState.failedOrbits = pState.failedOrbits + 1;
    end
end

