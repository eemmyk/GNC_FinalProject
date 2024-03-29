function [resultVector_o, paramVector_o, theta_super] = updateParameters(updateTOF, pSettings)
    global pState paramVector;

    if updateTOF
        global currentGuessVector
        %Initial guesses for parameters
        r1 = pSettings.a_initial;
        r2 = pSettings.a_final;

        currentGuessVector = [r1 r2 Inf 0, 0, 0, 0, 0];

        tofHandle = @(angle) fSolveTofFunction(angle, pSettings, pState);
        fminbnd(tofHandle, 0, 2.5*pi + pState.N*2*pi, pSettings.opt_tf_angle);

        r1 = currentGuessVector(7);
        r2 = currentGuessVector(8);
        bestTheta = currentGuessVector(4);
        paramVector.theta1 = currentGuessVector(5);
        paramVector.theta2 = currentGuessVector(6);

        TOF_estimation = (2*pState.N*pi + bestTheta)*sqrt(((r1+r2) * pSettings.TOF_corrMult)^3/(8*pSettings.mju)) + sqrt((abs(r1-r2))^3/(8*pSettings.mju));
        pState.tof_current = TOF_estimation;
        pState.initial_tof = TOF_estimation;
        paramVector.mju = pSettings.mju;

    end

    if pState.previousTime ~= pState.currentTime
        T_nu = mod(pState.currentTime - pSettings.Tp1, pSettings.P1);
        if T_nu ~= 0        
            n = pSettings.n1;
            e = pSettings.e1;
            nu_guess = paramVector.theta1 - pSettings.omega1;

            nu1_i = nuFromTime(T_nu, n, e, nu_guess);
        else
            nu1_i = 0;
        end
    
        theta1 = nu1_i + pSettings.omega1;
        gamma1 = asin(pSettings.e1 * sin(nu1_i) / sqrt(1+2*pSettings.e1*cos(nu1_i) + pSettings.e1^2));
    
        r1 = pSettings.p1 / (1+pSettings.e1*cos(nu1_i));
        theta1_dot = sqrt(pSettings.mju/pSettings.a_initial^3) * pSettings.a_initial^2/r1^2 * sqrt(1-pSettings.e1^2);
        
        T_nu = mod(pState.currentTime - pSettings.Tp2, pSettings.P2);
        if T_nu ~= 0
            n = pSettings.n2;
            e = pSettings.e2;
            nu_guess = paramVector.theta2 - pSettings.omega2;

            nu2_i = nuFromTime(T_nu, n, e, nu_guess);
        else
            nu2_i = 0;
        end
    else
        gamma1 = paramVector.gamma1;
        theta1_dot = paramVector.theta1_dot;
        r1 = paramVector.r1;
        theta1 = paramVector.theta1;
        nu2_i = paramVector.nu2_i;
    end
    pState.previousTime = pState.currentTime;
    
    T_nu = mod(pState.currentTime + pState.tof_current - pSettings.Tp2, pSettings.P2);
    if T_nu ~= 0
        n = pSettings.n2;
        e = pSettings.e2;
        nu_guess = paramVector.theta2 - pSettings.omega2;

        nu2 = nuFromTime(T_nu, n, e, nu_guess);
    else
        nu2 = 0;
    end

    theta2 = nu2 + pSettings.omega2;    
    theta_tilde = theta2 - theta1;

    if theta_tilde < 0
        theta_tilde = theta_tilde + 2*pi;
    elseif theta_tilde > 2*pi
        theta_tilde = theta_tilde - 2*pi;
    end
   

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
            geometricAngle = (geometricAngle + pi * (pSettings.safeTransferAngleMultiplier - 1)) / pSettings.safeTransferAngleMultiplier;
            if theta_f < geometricAngle
                theta_f = theta_f + 2*pi;
            end
        end
    end
    
    theta_vec = pSettings.approxSpace.*theta_f;
    
    %Pre power the long vector instead of doing it many times
    theta_vec1 = theta_vec;
    theta_vec2 = theta_vec1.*theta_vec;
    theta_vec3 = theta_vec2.*theta_vec;
    theta_vec4 = theta_vec3.*theta_vec;
    theta_vec5 = theta_vec4.*theta_vec;
    theta_vec6 = theta_vec5.*theta_vec;

    a = 1/r1;
    b = -tan(gamma1) / r1;
    c = 1/(2*r1) * (pSettings.mju / (r1^3 * theta1_dot^2) - 1);
    
    efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
                -48*theta_f     18*theta_f^2 -2*theta_f^3; 
                 20            -8*theta_f     theta_f^2];
    

    theta_super = [theta_vec1; theta_vec2; theta_vec3; theta_vec4; theta_vec5; theta_vec6];

    paramVector.gamma1 = gamma1;
    paramVector.gamma2 = gamma2;
    paramVector.theta_f = theta_f;
    paramVector.theta1_dot = theta1_dot;
    paramVector.theta2_dot = theta2_dot;
    paramVector.r1 = r1;
    paramVector.r2 = r2;
    paramVector.theta1 = theta1;
    paramVector.theta2 = theta2;
    paramVector.nu2_i = nu2_i;
    paramVector.r2_i = r2_i;
    paramVector.a = a;
    paramVector.b = b;
    paramVector.c = c;
    paramVector.efg_Mat_1 = efg_Mat_1;

    %% Check if TOF is a solution 
    

    %Limit the distances geometrically
    if theta_f-gamma2+gamma1 < pi
        %Remember the multiplier for theta_f is a test
        objectDistance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta_f));
      
        alpha1 = acos((2*r1^2 - 2*r1*r2*cos(theta_f)) / (2*r1*objectDistance));

        coord_A = tan(theta1-gamma1+pi/2);
        coord_B = tan(theta2+gamma2+pi/2);
        coord_C = -tan(theta1-gamma1+pi/2)*cos(theta1)*r1+sin(theta1)*r1;
        coord_D = -tan(theta2+gamma2+pi/2)*cos(theta2)*r2+sin(theta2)*r2;
        
        x_max = (coord_D - coord_C) / (coord_A - coord_B);
        y_max = coord_A * (coord_D - coord_C) / (coord_A - coord_B) + coord_C;

        altGeometricMaxRadius = sqrt(x_max^2 + y_max^2);

        geometricMaxRadius = min(altGeometricMaxRadius, pSettings.rMax);

        geometricMinRadius = sin(alpha1) * r1;
        geometricMinRadius = max(geometricMinRadius, pSettings.rMin);


%         figure(3)
%         hold on
%         rectangle('Position',[-geometricMinRadius, -geometricMinRadius, 2*geometricMinRadius, 2*geometricMinRadius],'Curvature',[1 1])
%         rectangle('Position',[-geometricMaxRadius, -geometricMaxRadius, 2*geometricMaxRadius, 2*geometricMaxRadius],'Curvature',[1 1])


        
    else
        geometricMaxRadius = pSettings.rMax;
        geometricMinRadius = pSettings.rMin;
    end

    % Could implement real root finding algorithm for this part.



    %Calculate the minimum value for d
    d_minimum = fFindRadiusFunction(paramVector, geometricMaxRadius);

    %endRanges = [theta_super(:,1:ceil(pSettings.intApprox * 0.2)), theta_super(:,floor(pSettings.intApprox*0.8):end)];
    
    [timePart, radiusPart] = fTimeMinReal(theta_super, d_minimum, paramVector, pSettings.a_initial);

    if (timePart < 0) || (radiusPart < 0)
        y0_1 = timePart;
        y0_2 = radiusPart;
        x0 = d_minimum;
        
        if d_minimum < 0
            d_minimum = d_minimum / pSettings.dAdjustment + 1e-15;
        else
            d_minimum = d_minimum * pSettings.dAdjustment + 1e-15;
        end

        [timePart, radiusPart] = fTimeMinReal(theta_super, d_minimum, paramVector, pSettings.a_initial);

        if (timePart < 0) || (radiusPart < 0)
            y1_1 = timePart;
            y1_2 = radiusPart;
            x1 = d_minimum;
    
            k1 = (y1_1-y0_1)/(x1-x0);
            k2 = (y1_2-y0_2)/(x1-x0);
    
            %Where does the troublesome time part reach 0.1
            d_minimum_1 = (0.1 - y1_1)/k1 + x1;
            d_minimum_2 = (0.1 - y1_2)/k2 + x1;
    
            if radiusPart < 0
                d_minimum = max(d_minimum_1, d_minimum_2);
            else
                d_minimum = d_minimum_1;
            end
        end
    end

    %Calculate the maximum value for d
    d_maximum = fFindRadiusFunction(paramVector, geometricMinRadius);

    %middleRange = theta_super(:,floor(pSettings.intApprox*0.4):ceil(pSettings.intApprox*0.6));

    [timePart, radiusPart] = fTimeMinReal(theta_super, d_maximum, paramVector, pSettings.a_initial);

    if (timePart < 0) || (radiusPart < 0)
        y0_1 = timePart;
        %y0_2 = radiusPart;
        x0 = d_maximum;
        
        if d_maximum < 0
            d_maximum = d_maximum * pSettings.dAdjustment - 1e-15;
        else
            d_maximum = d_maximum / pSettings.dAdjustment - 1e-15;
        end

        [timePart, ~] = fTimeMinReal(theta_super, d_maximum, paramVector, pSettings.a_initial);

        if (timePart < 0) || (radiusPart < 0)
            y1_1 = timePart;
            x1 = d_maximum;
    
            k1 = (y1_1-y0_1)/(x1-x0);
            %k2 = (y1_2-y0_2)/(x1-x0);
    
            %Where does the troublesome time part reach 0.1
            d_maximum_1 = (0.1 - y1_1)/k1 + x1;
    
            d_maximum = d_maximum_1;
        end
    end

    [timePart, radiusPart] = fTimeMinReal(theta_super, (d_minimum + d_maximum) / 2, paramVector, pSettings.a_initial);

    realOrbit = 0;
    if (timePart > 0) && (radiusPart > 0)
        realOrbit = 1;
    end

    resultVector = [d_minimum, d_maximum, realOrbit];
    resultVector_o = resultVector;
    paramVector_o = paramVector;
end

%% Function for finding the d-coefficient where orbit reaches r (at theta_f/2)
function [d_sol] = fFindRadiusFunction(paramVector, rTarget)

    mju = paramVector.mju;
    gamma1 = paramVector.gamma1;
    gamma2 = paramVector.gamma2;
    theta_f = paramVector.theta_f;
    theta1_dot = paramVector.theta1_dot;
    theta2_dot = paramVector.theta2_dot;
    r1 = paramVector.r1;
    r2 = paramVector.r2;

    d_sol = (64*(1/rTarget - 1/r1 + (theta_f^6*(((5*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (10*tan(gamma1)*theta_f)/r1 + 10/r1 - 10/r2)/theta_f^6 - ((4*tan(gamma2))/r2 - (4*tan(gamma1))/r1 + (4*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^5 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^4)))/64 + (theta_f^4*(((15*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/(2*r1) - (15*tan(gamma1)*theta_f)/r1 + 15/r1 - 15/r2)/theta_f^4 - ((5*tan(gamma2))/r2 - (5*tan(gamma1))/r1 + (5*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^3 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^2)))/16 - (theta_f^5*(((12*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (24*tan(gamma1)*theta_f)/r1 + 24/r1 - 24/r2)/theta_f^5 - ((9*tan(gamma2))/r2 - (9*tan(gamma1))/r1 + (9*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^4 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/theta_f^3))/32 + (theta_f*tan(gamma1))/(2*r1) - (theta_f^2*(mju/(r1^3*theta1_dot^2) - 1))/(8*r1)))/theta_f^3;
 
end

%% Solve nu around an orbit at a given time T
function [angleError] = nuSolver(nu_time, T, n, e)
    if nu_time <= pi
        E = 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    else
        E = 2*pi + 2*atan(tan(nu_time/2)/sqrt((1+e)/(1-e)));
    end

    angleError = E-e*sin(E) - n*T;
end

%% Solve nu around an orbit at a given time T
function [nu] = nuFromTime(T, n, e, nu_guess)

    M = n*T;
    nu = M;

    if e ~= 0
    
        root = sqrt(-e^3 * (-8 + 24*e - 24*e^2 + 8*e^3 - 9*e*M^2));
    
        nom = -2*e + 2*e^2 + (3*M*e^2 + root)^(2/3);
        denom = e * (3*M*e^2 + root)^(1/3);
    
        E =  nom/denom;
        
        if E <= pi
            nu = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
        else
            nu = 2*pi + 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
        end
    
        if M <= pi
            E2 = 2*atan(tan(M/2)/sqrt((1+e)/(1-e)));
        else
            E2 = 2*pi + 2*atan(tan(M/2)/sqrt((1+e)/(1-e)));
        end


        xnm1 = M;
        fnm1 = E2 - e*sin(E2) - M;
        xn = nu;
        fn = E - e*sin(E) - M;


       
        iter = 0;
        while (abs(fn) > 0.001) && (iter < 100)
            xnp1 = mod(xn - fn * (xn-xnm1) / (fn-fnm1), 2*pi);
            fnm1 = fn;
            xnm1 = xn;
            xn = xnp1;
            fn = nuSolver(xnp1, T, n, e);

            iter = iter + 1;

            if iter == 50
                if abs(fnm1) < abs(fn)
                    fn = nuSolver(nu_guess, T, n, e);
                    xn = nu_guess;
                else
                    fnm1 = nuSolver(nu_guess, T, n, e);
                    xnm1 = nu_guess;            
                end
            end
%             if abs(fn) < abs(minFn)
%                 minFn = fn;
%                 bestXn = xn;
%             end
        end

        nu = xn;%bestXn;

    end

end


%% Solve the bounds of the d-coefficient from part of the time function
function [timeImgPart, radiusImgPart] = fTimeMinReal(theta, d, paramVector, a_initial)

    theta_f = paramVector.theta_f;
    r2 = paramVector.r2;
    
    a = paramVector.a;
    b = paramVector.b;
    c = paramVector.c;
                
    efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
                -tan(paramVector.gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
                paramVector.mju/(r2^4*paramVector.theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
    
    efg = 1/(2*theta_f^6) * paramVector.efg_Mat_1 * efg_Mat_2;
    
    e = efg(1);
    f = efg(2);
    g = efg(3);

    radiusImgPart = min(a_initial*(a + b.*theta(1,:) + c.*theta(2,:) + d.*theta(3,:) + e.*theta(4,:) + f.*theta(5,:) + g.*theta(6,:)));
    timeImgPart = min(a_initial*(a + 2.*c + (6.*d + b).*theta(1,:) + (12.*e + c).*theta(2,:) + (20.*f + d).*theta(3,:) + (30.*g + e).*theta(4,:) + f.*theta(5,:) + g.*theta(6,:)));
end

%% For Solving the initial TOF guess 
function [meetAngleError] = fSolveTofFunction(transferAngle, pSettings, pState)
    global currentGuessVector

    r1 = currentGuessVector(1);
    r2 = currentGuessVector(2);
    bestAngleError = currentGuessVector(3);
   
    TOF_estimation = (2*pState.N*pi + transferAngle)*sqrt(((r1+r2) * pSettings.TOF_corrMult)^3/(8*pSettings.mju)) + sqrt((abs(r1-r2))^3/(8*pSettings.mju));
   
    T_nu = mod(pState.currentTime - pSettings.Tp1, pSettings.P1);
    if T_nu ~= 0
        nuHandle = @(angle) nuSolver(angle, T_nu, pSettings.n1, pSettings.e1);

        nu1_i = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
    else
        nu1_i = 0;
    end

    theta1 = nu1_i + pSettings.omega1;
    gamma1 = asin(pSettings.e1 * sin(nu1_i) / sqrt(1+2*pSettings.e1*cos(nu1_i) + pSettings.e1^2));

    r1 = pSettings.p1 / (1+pSettings.e1*cos(nu1_i));

    T_nu = mod(pState.currentTime + TOF_estimation - pSettings.Tp2, pSettings.P2);
    if T_nu ~= 0
        nuHandle = @(angle) nuSolver(angle, T_nu, pSettings.n2, pSettings.e2);

        nu2 = fzero(nuHandle, [0, 2*pi], pSettings.opt_nu_fzero);
    else
        nu2 = 0;
    end

    r2 = pSettings.p2 / (1+pSettings.e2*cos(nu2));    
    gamma2 = asin(pSettings.e2 * sin(nu2) / sqrt(1+2*pSettings.e2*cos(nu2) + pSettings.e2^2));

    theta2 = nu2 + pSettings.omega2;
    theta_tilde = theta2 - theta1;

    if theta_tilde < 0
        theta_tilde = theta_tilde + 2*pi;
    elseif theta_tilde > 2*pi
        theta_tilde = theta_tilde - 2*pi;
    end
    
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
        geometricAngle = (geometricAngle + pi * (pSettings.safeTransferAngleMultiplier - 1)) / pSettings.safeTransferAngleMultiplier;
        if theta_f < geometricAngle
            theta_f = theta_f + 2*pi;
        end
    end

    meetAngleError = abs(transferAngle - theta_f);

    if meetAngleError < bestAngleError
        currentGuessVector(3) = meetAngleError;
        currentGuessVector(4) = transferAngle;
        currentGuessVector(5) = theta1;
        currentGuessVector(6) = theta2;
        currentGuessVector(7) = r1;
        currentGuessVector(8) = r2;
    end


    currentGuessVector(1) = r1;
    currentGuessVector(2) = r2;

end
