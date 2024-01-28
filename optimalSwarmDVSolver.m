function [resultMatrix, dateResultMatrix, tofResultMatrix, paramResultMatrix, efgResultMatrix] = optimalSwarmDVSolver(swarmParams, pSettings, twMapInds)
    global pState deltaResult

    leadNuPointCountPerOrbit = swarmParams.leadNuPointCountPerOrbit;
    datePointCount = swarmParams.datePointCount;
    extraOrbitChecks = swarmParams.extraOrbitChecks;
    %deployNuVector = swarmParams.deployNuVector;
    %targetNuVector = swarmParams.targetNuVector;
    nuCount = swarmParams.nuCount;
    dateVector = swarmParams.dateVector;
    %leadNuVector = swarmParams.leadNuVector;
    tofMatrix = swarmParams.tofMatrix;
    deployTime = swarmParams.deployTime;
    targetTime = swarmParams.targetTime;
    perigeeTimeVectorInitial = swarmParams.perigeeTimeVectorInitial;
    perigeeTimeVectorFinal = swarmParams.perigeeTimeVectorFinal;

    resultMatrix = zeros(nuCount, nuCount);
    dateResultMatrix = zeros(nuCount, nuCount);
    tofResultMatrix = zeros(nuCount, nuCount);
    %leadNuResultMatrix = zeros(nuCount, nuCount);
    paramResultMatrix = zeros(nuCount, nuCount, 16);
    efgResultMatrix = zeros(nuCount, nuCount, 9);

    for swarmIndTarget = 1:nuCount
        %nuEnd = targetNuVector(swarmIndTarget);
        for swarmIndDeploy = 1:nuCount
            fprintf("Finding orbits: <strong>%.1f%%</strong>\n", 100 * ((swarmIndTarget-1) * nuCount + swarmIndDeploy)/(nuCount * nuCount));
            %nuStart = deployNuVector(swarmIndDeploy);
            
            Tp1 = perigeeTimeVectorInitial(swarmIndDeploy);
            Tp2 = perigeeTimeVectorFinal(swarmIndTarget);

            deltaResult = Inf;

            for dateInd = 1:datePointCount
                date = dateVector(dateInd);
                for tofInd = 1:leadNuPointCountPerOrbit * (extraOrbitChecks+1)
                    tof = tofMatrix(swarmIndTarget, dateInd, tofInd);
                    
%                     leadNuIndex = mod(tofInd, leadNuPointCountPerOrbit);
%                     if leadNuIndex == 0
%                         leadNuIndex = leadNuPointCountPerOrbit;
%                     end

                    %leadNu = leadNuVector(leadNuIndex);


                    pState.previousTime = -1;

                    if pSettings.solveDate == 1
                        pState.currentTime = date;
                        pState.tof_current = tof;
                    else
                        pState.tof_current = tof;
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
                
                    
                    d_solution = 0;
                    deltaV_o = 1e24; %A big number
                
                    for N_current = N_start:N_end
                        
                        pState.testedOrbits = pState.testedOrbits + 1;
                        
                        pState.N = N_current;

                        if swarmIndDeploy == 1
                            0;
                        end

                        [resultVector, paramVector, theta_super] = updateSwarmParameters(Tp1, Tp2, deployTime, targetTime, 0, pSettings);
                        dT = theta_super(1,2) - theta_super(1,1);
                        tof_current = pState.tof_current;
                
                        d_minimum = resultVector(1);
                        d_maximum = resultVector(2);
                        realOrbit = resultVector(3);
                
                         if ~realOrbit
                            continue
                         end
                    
                        %Calculate minimum and maximum time of flight
                        x_min = d_minimum;
                        f_min = fTimeFunction(x_min, theta_super, dT, paramVector) - tof_current;
                        x_max = d_maximum;
                        f_max = fTimeFunction(x_max, theta_super, dT, paramVector) - tof_current;
                                        
                        crossing = (f_min < 0) ~= (f_max < 0);
                
                        if ~crossing %|| imaginary
                            continue
                        end
                        
                        d_minimum_in = d_minimum;
                        d_maximum_in = d_maximum;
                        
                        if (d_minimum < 0) ~= (d_maximum < 0)
                            f_zero = fTimeFunction(0, theta_super, dT, paramVector) - tof_current;
                %             if ~isreal(f_zero)
                %                 continue
                %             else
                            if (f_zero > 0) && (f_min > 0)
                                d_minimum_in = 0;
                                f_min = f_zero;
                            else
                                d_maximum_in = 0;
                                f_max = f_zero;
                            end
                        end
                
                        crossing = (f_min < 0) ~= (f_max < 0);
                        if ~crossing
                            continue
                        end
                
                
                        % Customisted to be specific to transferTimeSolution-function
                        % Parameters are passed directly to the function instead of through a handle
                        d_solution = fMyFastZero(@transferTimeSolution, [d_minimum_in, d_maximum_in], [f_min, f_max], pSettings.opt_tof_fzero, paramVector, pState.tof_current, theta_super);
                
                        deltaV_o = fJerkFunction(d_solution, theta_super, dT, paramVector);
                        trueSolution = 1;
                        pState.successfulOrbits = pState.successfulOrbits + 1;
                
                        if deltaV_o <= localBestDV
                            localBestDV = deltaV_o;
                        else
                            break
                        end 
                
                        if localBestDV < deltaResult
                            %Save best results as globals
                            global dateOptimal tof_optimal d_opt;
                            global paramVector_opt resultVector_opt
                    
                            deltaResult = localBestDV;
                            d_opt = d_solution;
                    
                            dateOptimal = pState.currentTime;
                            tof_optimal = pState.tof_current;

                            %leadNuOptimal = leadNu;
                    
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
                    
                end
            end

            resultMatrix(swarmIndDeploy, swarmIndTarget) = deltaResult;
            dateResultMatrix(swarmIndDeploy, swarmIndTarget) = dateOptimal;
            tofResultMatrix(swarmIndDeploy, swarmIndTarget) = tof_optimal;
            efgResultMatrix(swarmIndDeploy, swarmIndTarget, :) = [paramVector_opt.efg_Mat_1(1, :), paramVector_opt.efg_Mat_1(2, :), paramVector_opt.efg_Mat_1(3, :)];
            %leadNuResultMatrix(swarmIndDeploy, swarmIndTarget) = leadNuOptimal;

            mju_result = paramVector_opt.mju;
            gamma1_result = paramVector_opt.gamma1;
            gamma2_result = paramVector_opt.gamma2;
            theta_f_result = paramVector_opt.theta_f;
            theta1_dot_result = paramVector_opt.theta1_dot;
            theta2_dot_result = paramVector_opt.theta2_dot;
            r1_result = paramVector_opt.r1;
            r2_result = paramVector_opt.r2;
            theta1_result = paramVector_opt.theta1;
            theta2_result = paramVector_opt.theta2;
            nu2_i_result = paramVector_opt.nu2_i;
            r2_i_result = paramVector_opt.r2_i;
            a_result = paramVector_opt.a;
            b_result = paramVector_opt.b;
            c_result = paramVector_opt.c;
            d_result = d_opt;

            paramResultMatrix(swarmIndDeploy, swarmIndTarget, :) = [mju_result, gamma1_result, gamma2_result, theta_f_result, theta1_dot_result, theta2_dot_result, r1_result, r2_result, theta1_result, theta2_result, nu2_i_result, r2_i_result, a_result, b_result, c_result, d_result];
   
        end
    end
end