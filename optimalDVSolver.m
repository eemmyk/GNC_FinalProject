function [deltaV_o] = optimalDVSolver(inputVec, pSettings)
    global d_solution paramVector resultVector deltaResult theta_vec pState; 

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



% function [timeStep] = fTimeDervivative(d, theta, paramVector)
% 
%     mju = paramVector(1);
%     gamma1 = paramVector(2);
%     gamma2 = paramVector(3);
%     theta_f = paramVector(4);
%     theta1_dot = paramVector(5);
%     theta2_dot = paramVector(6);
%     r1 = paramVector(7);
%     r2 = paramVector(8);
% 
% %     a = 1/r1;
% %     b = -tan(gamma1) / r1;
% %     c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);
% %     
% %     efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
% %                 -48*theta_f     18*theta_f^2 -2*theta_f^3; 
% %                  20            -8*theta_f     theta_f^2];
% %     
% %     efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
% %                 -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
% %                 mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
% %     
% %     efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
% %     
% %     e = efg(1);
% %     f = efg(2);
% %     g = efg(3);
%      
%     timeStep = ((6.*d - tan(gamma1)./r1 - 2.*theta.*((6.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^2 + (180.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (60.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) + 3.*d.*theta.^2 - 4.*theta.^3.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) + (15.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (5.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) - 6.*theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) + (10.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (4.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + 5.*theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 + (24.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (9.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) - 4.*theta.^3.*((15.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^4 + (300.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (120.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + 3.*theta.^2.*((20.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^3 + (480.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (180.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) + (theta.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./(mju.*(theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) - theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) - theta.^6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) + d.*theta.^3 + 1./r1 - (theta.*tan(gamma1))./r1 + (theta.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)).^4) + (4.*(tan(gamma1)./r1 - 3.*d.*theta.^2 + 4.*theta.^3.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) + (15.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (5.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) + 6.*theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) + (10.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (4.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) - 5.*theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 + (24.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (9.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) - (theta.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1).*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta + d.*theta.^3 + 1./r1 - theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) + (15.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (5.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) - theta.^6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) + (10.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (4.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 + (24.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (9.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) - theta.^2.*((6.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^2 + (180.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (60.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) - theta.^4.*((15.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^4 + (300.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (120.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + theta.^3.*((20.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^3 + (480.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (180.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) - (theta.*tan(gamma1))./r1 + (theta.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./(mju.*(theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) - theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) - theta.^6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) + d.*theta.^3 + 1./r1 - (theta.*tan(gamma1))./r1 + (theta.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)).^5))./(2.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta - theta.^6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) - theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) + theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) + d.*theta.^3 - theta.^2.*(((6.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 36.*d.*theta_f + 6./r2 - (6.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^2 - (180.*d.*theta_f.^2 + (60.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (60.*tan(gamma1))./r1 + (60.*tan(gamma2))./r2)./theta_f.^3 + (180.*d.*theta_f.^3 + (90.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (180.*tan(gamma1).*theta_f)./r1 + 180./r1 - 180./r2)./theta_f.^4) - theta.^4.*(((15.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 90.*d.*theta_f + 15./r2 - (15.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^4 - (360.*d.*theta_f.^2 + (120.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (120.*tan(gamma1))./r1 + (120.*tan(gamma2))./r2)./theta_f.^5 + (300.*d.*theta_f.^3 + (150.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (300.*tan(gamma1).*theta_f)./r1 + 300./r1 - 300./r2)./theta_f.^6) + theta.^3.*(((20.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 120.*d.*theta_f + 20./r2 - (20.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (540.*d.*theta_f.^2 + (180.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (180.*tan(gamma1))./r1 + (180.*tan(gamma2))./r2)./theta_f.^4 + (480.*d.*theta_f.^3 + (240.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (480.*tan(gamma1).*theta_f)./r1 + 480./r1 - 480./r2)./theta_f.^5) + 1./r1 - (theta.*tan(gamma1))./r1 + (theta.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1))./(mju.*(theta.^5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - ((9.*tan(gamma2))./r2 - (9.*tan(gamma1))./r1 + 27.*d.*theta_f.^2 + (9.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^4 + (24.*d.*theta_f.^3 + 24./r1 - 24./r2 - (24.*theta_f.*tan(gamma1))./r1 + (12.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^5) - theta.^4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - ((5.*tan(gamma2))./r2 - (5.*tan(gamma1))./r1 + 15.*d.*theta_f.^2 + (5.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^3 + (15.*d.*theta_f.^3 + 15./r1 - 15./r2 - (15.*theta_f.*tan(gamma1))./r1 + (15.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1))./theta_f.^4) - theta.^6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - ((4.*tan(gamma2))./r2 - (4.*tan(gamma1))./r1 + 12.*d.*theta_f.^2 + (4.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^5 + (10.*d.*theta_f.^3 + 10./r1 - 10./r2 - (10.*theta_f.*tan(gamma1))./r1 + (5.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^6) + d.*theta.^3 + 1./r1 - (theta.*tan(gamma1))./r1 + (theta.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)).^4)).^(1./2));
%  
% end
