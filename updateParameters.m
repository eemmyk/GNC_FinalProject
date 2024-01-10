function [] = updateParameters(updateTOF, pSettings)
    global theta_vec TOF_estimation

    global pState;

    %Outgoing globals
    global resultVector paramVector

    %Internally outgoing
    global r1 r2 theta1 theta2 theta_f gamma1 gamma2 bestAngleError bestTheta

    if updateTOF
        %Initial guesses for parameters
        r1 = pSettings.a_initial;
        r2 = pSettings.a_final;
        gamma1 = 0;
        gamma2 = 0;

        bestAngleError = Inf;

        tofHandle = @(angle) fSolveTofFunction(angle, pSettings, pState);

        fminbnd(tofHandle, 0, 2.5*pi, pSettings.opt_tf_angle);

        TOF_estimation = (2*pState.N*pi + bestTheta)*sqrt(((r1+r2) * pSettings.TOF_corrMult)^3/(8*pSettings.mju)) + sqrt((abs(r1-r2))^3/(8*pSettings.mju));
        pState.tof_current = TOF_estimation;
    end

    if pState.previousTime ~= pState.currentTime
        T_nu = mod(pState.currentTime - pSettings.Tp1, pSettings.P1);
        if T_nu ~= 0        
            nuHandle = @(angle) nuFromTime(angle, T_nu, pSettings.n1, pSettings.e1);

            nu1_i = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
        else
            nu1_i = 0;
        end
    
        theta1 = nu1_i + pSettings.omega1;
        gamma1 = asin(pSettings.e1 * sin(nu1_i) / sqrt(1+2*pSettings.e1*cos(nu1_i) + pSettings.e1^2));
    
        r1 = pSettings.p1 / (1+pSettings.e1*cos(nu1_i));
        theta1_dot = sqrt(pSettings.mju/pSettings.a_initial^3) * pSettings.a_initial^2/r1^2 * sqrt(1-pSettings.e1^2);
        
        T_nu = mod(pState.currentTime - pSettings.Tp2, pSettings.P2);
        if T_nu ~= 0
            nuHandle = @(angle) nuFromTime(angle, T_nu, pSettings.n2, pSettings.e2);
        
            nu2_i = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
        else
            nu2_i = 0;
        end
    else
        gamma1 = paramVector(2);
        theta1_dot = paramVector(5);
        r1 = paramVector(7);
        theta1 = paramVector(9);
        nu2_i = paramVector(11);

    end
    pState.previousTime = pState.currentTime;
    
    T_nu = mod(pState.currentTime + pState.tof_current - pSettings.Tp2, pSettings.P2);
    if T_nu ~= 0
        nuHandle = @(angle) nuFromTime(angle, T_nu, pSettings.n2, pSettings.e2);

        nu2 = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
    else
        nu2 = 0;
    end

    theta2 = nu2 + pSettings.omega2;    
    theta_tilde = theta2 - theta1;

    gamma2 = asin(pSettings.e2 * sin(nu2) / sqrt(1+2*pSettings.e2*cos(nu2) + pSettings.e2^2));
    r2 = pSettings.p2 / (1+pSettings.e2*cos(nu2));    
    r2_i = pSettings.p2 / (1+pSettings.e2*cos(nu2_i));    
    theta2_dot = sqrt(pSettings.mju./pSettings.a_final.^3).*pSettings.a_final.^2./r2.^2.*sqrt(1-pSettings.e2.^2);

    theta_f = 2.*pi.*pState.N + theta_tilde;
    
    %This is (less) fucked
    if theta_f < pi
        cutDist1 = r1;
        cutDist2 = r2;
        if gamma1 < 0
            cutDist2 = sin(pi/2 + gamma1)*r1/sin(pi/2-gamma1-theta_f);
        end
    
        if gamma2 > 0
            cutDist1 = sin(pi/2 - gamma2)*r2/sin(pi/2+gamma2-theta_f);
        end
        
        if ((0 < cutDist1) && (cutDist1 < r1)) || ((0 < cutDist2) && (cutDist2 < r2))
            theta_f = theta_f + 2*pi;
        else
            geometricAngle = pi/2 - asin(min(r1,r2)/max(r1,r2));
            geometricAngle = (2 * geometricAngle - pi) / pSettings.safeTransferAngleMultiplier + pi;
            if theta_f < geometricAngle * pSettings.safeTransferAngleMultiplier
                theta_f = theta_f + 2*pi;
            end
        end
    end
    
    theta_vec = linspace(pSettings.theta_0, theta_f, pSettings.intApprox);
    theta_vec_acc = linspace(pSettings.theta_0, theta_f, pSettings.intApprox);

    paramVector = [pSettings.mju, gamma1, gamma2, theta_f, theta1_dot, theta2_dot, r1, r2, theta1, theta2, nu2_i, r2_i];

    %% Check if TOF is a solution 

    %Limit the distances geometrically
    if theta_f < pi
        objectDistance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta_f));
      
        alpha1 = acos((2*r1^2 - 2*r1*r2*cos(theta_f)) / (2*r1*objectDistance));

        objectToMax2 = sin((pi/2)-alpha1) * objectDistance / sin(2*pi - theta_f);

        geometricMaxRadius = sqrt(objectToMax2^2 + r2^2);
        geometricMaxRadius = min(geometricMaxRadius, pSettings.rMax);

        geometricMinRadius = sin(alpha1) * r1;
        geometricMinRadius = max(geometricMinRadius, pSettings.rMin);
    else
        geometricMaxRadius = pSettings.rMax;
        geometricMinRadius = pSettings.rMin;
    end

    %Calculate the minimum value for d
    d_minimum = fFindRadiusFunction(paramVector, geometricMaxRadius);

    endRanges = [theta_vec_acc(floor(pSettings.intApprox*0.05):ceil(pSettings.intApprox * 0.1)), theta_vec_acc(floor(pSettings.intApprox*0.9):ceil(pSettings.intApprox * 0.95))];

    tof_max = min(fTimeMinReal(endRanges, d_minimum, paramVector, pSettings.a_initial));
    while tof_max < 0
        %Move a relative amount towards larger values and an arbitary
        %small step to cross zero properly
        if d_minimum < 0
            d_minimum = d_minimum / pSettings.dAdjustment + 1e-15;
        else
            d_minimum = d_minimum * pSettings.dAdjustment + 1e-15;
        end

        tof_max = min(fTimeMinReal(endRanges, d_minimum, paramVector, pSettings.a_initial));
        %pState.No real solutions exist for orbit
        if d_minimum > 1
            d_maximum = 1;
            return
        end
    end
    
    %Calculate the maximum value for d
    d_maximum = fFindRadiusFunction(paramVector, geometricMinRadius);

    middleRange = theta_vec_acc(floor(pSettings.intApprox*0.45):ceil(pSettings.intApprox*0.55));

    tof_min = min(fTimeMinReal(middleRange, d_maximum, paramVector, pSettings.a_initial));
    while tof_min < 0
        %Move a relative amount towards larger values and an arbitary
        %small step to cross zero properly
        if d_maximum < 0
            d_maximum = d_maximum * pSettings.dAdjustment - 1e-15;
        else
            d_maximum = d_maximum / pSettings.dAdjustment - 1e-15;
        end

        tof_min = min(fTimeMinReal(middleRange, d_maximum, paramVector, pSettings.a_initial));
    end

    resultVector = [d_minimum, d_maximum, TOF_estimation];
end

%% Function for finding the d-coefficient where orbit reaches minimum r
function [d_sol] = fFindRadiusFunction(paramVector, rTarget)

    mju = paramVector(1);
    gamma1 = paramVector(2);
    gamma2 = paramVector(3);
    theta_f = paramVector(4);
    theta1_dot = paramVector(5);
    theta2_dot = paramVector(6);
    r1 = paramVector(7);
    r2 = paramVector(8);

    d_sol = (64*(1/rTarget - 1/r1 + (theta_f^6*(((5*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (10*tan(gamma1)*theta_f)/r1 + 10/r1 - 10/r2)/theta_f^6 - ((4*tan(gamma2))/r2 - (4*tan(gamma1))/r1 + (4*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^5 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^4)))/64 + (theta_f^4*(((15*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/(2*r1) - (15*tan(gamma1)*theta_f)/r1 + 15/r1 - 15/r2)/theta_f^4 - ((5*tan(gamma2))/r2 - (5*tan(gamma1))/r1 + (5*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^3 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^2)))/16 - (theta_f^5*(((12*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (24*tan(gamma1)*theta_f)/r1 + 24/r1 - 24/r2)/theta_f^5 - ((9*tan(gamma2))/r2 - (9*tan(gamma1))/r1 + (9*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^4 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/theta_f^3))/32 + (theta_f*tan(gamma1))/(2*r1) - (theta_f^2*(mju/(r1^3*theta1_dot^2) - 1))/(8*r1)))/theta_f^3;
 
end

%% Solve nu around an orbit at a given time T
function [angleError] = nuFromTime(nu_time, T, n, e)
    if nu_time <= pi
        E = 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    else
        E = 2*pi + 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    end

    angleError = E-e*sin(E) - n*T;
end

%% Solve maximum d-coefficient from part of the time function
function [timeImgPart] = fTimeMinReal(theta, d, paramVector, a_initial)
    
    mju = paramVector(1);
    gamma1 = paramVector(2);
    gamma2 = paramVector(3);
    theta_f = paramVector(4);
    theta1_dot = paramVector(5);
    theta2_dot = paramVector(6);
    r1 = paramVector(7);
    r2 = paramVector(8);

    
    a = 1/r1;
    b = -tan(gamma1) / r1;
    c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);
    
    efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
                -48*theta_f     18*theta_f^2 -2*theta_f^3; 
                 20            -8*theta_f     theta_f^2];
    
    efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
                -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
                mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
    
    efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
    
    e = efg(1);
    f = efg(2);
    g = efg(3);
    
    r = 1 ./ (a + b.*theta + c.*theta.^2 + d.*theta.^3 + e.*theta.^4 + f.*theta.^5 + g.*theta.^6);

    timeImgPart = a_initial*(1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4);
end

%% For Solving the initial TOF guess 
function [meetAngleError] = fSolveTofFunction(transferAngle, pSettings, pState)
    global TOF_estimation r1 r2 theta1 theta2 theta_f gamma1 gamma2 bestAngleError bestTheta

%     geometricAngle = pi/2 - asin(min(r1,r2)/max(r1,r2));
%     if transferAngle < geometricAngle * pSettings.safeTransferAngleMultiplier
%         transferAngle = transferAngle + 2*pi;
%     end

    TOF_estimation = (2*pState.N*pi + transferAngle)*sqrt(((r1+r2) * pSettings.TOF_corrMult)^3/(8*pSettings.mju)) + sqrt((abs(r1-r2))^3/(8*pSettings.mju));
   
    T_nu = mod(pState.currentTime - pSettings.Tp1, pSettings.P1);
    if T_nu ~= 0
        nuHandle = @(angle) nuFromTime(angle, T_nu, pSettings.n1, pSettings.e1);

        nu1_i = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
    else
        nu1_i = 0;
    end

    theta1 = nu1_i + pSettings.omega1;
    gamma1 = asin(pSettings.e1 * sin(nu1_i) / sqrt(1+2*pSettings.e1*cos(nu1_i) + pSettings.e1^2));

    r1 = pSettings.p1 / (1+pSettings.e1*cos(nu1_i));

    T_nu = mod(pState.currentTime + TOF_estimation - pSettings.Tp2, pSettings.P2);
    if T_nu ~= 0
        nuHandle = @(angle) nuFromTime(angle, T_nu, pSettings.n2, pSettings.e2);

        nu2 = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
    else
        nu2 = 0;
    end
    
    gamma2 = asin(pSettings.e2 * sin(nu2) / sqrt(1+2*pSettings.e2*cos(nu2) + pSettings.e2^2));

    theta2 = nu2 + pSettings.omega2;
    theta_tilde = theta2 - theta1;

    r2 = pSettings.p2 / (1+pSettings.e2*cos(nu2));    
    
    theta_f = 2.*pi.*pState.N + theta_tilde;

    cutDist1 = r1;
    cutDist2 = r2;
    if gamma1 < 0
        cutDist2 = sin(pi/2 + gamma1)*r1/sin(pi/2-gamma1-theta_f);
    end

    if gamma2 > 0
        cutDist1 = sin(pi/2 - gamma2)*r2/sin(pi/2+gamma2-theta_f);
    end

    if ((0 < cutDist1) && (cutDist1 < r1)) || ((0 < cutDist2) && (cutDist2 < r2))
        theta_f = theta_f + 2*pi;
    else
        geometricAngle = pi/2 - asin(min(r1,r2)/max(r1,r2));
        geometricAngle = (2 * geometricAngle - pi) / pSettings.safeTransferAngleMultiplier + pi;
        if theta_f < geometricAngle
            theta_f = theta_f + 2*pi;
        end
    end
    
    meetAngleError = abs(transferAngle - theta_f);

    if meetAngleError < bestAngleError
        bestAngleError = meetAngleError;
        bestTheta = theta_f;
    end
end
