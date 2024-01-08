function [] = updateParameters(updateTOF)
    % Calculating the orbital parameters
    %initial orbit parameters
    global a_initial a_final currentTime mju N tof_current;
    global omega1 omega2 e1 e2;
    global theta_vec;
    global intApprox plotAccuracy;
    global Tp1 Tp2 TOF_estimation;
    global d_minimum d_maximum rMin rMax;
    global theta_0;
    global n1 P1 n2 P2 p1 p2;
    global previousTime unfeasibleOrbit safeTransferAngle;
    global TOF_corrMult dAdjustment;

    global paramVector;

    global theta_vec_acc;
    
    global opt_nu_fzero opt_tf_angle opt_d_lim_fzero;

    if updateTOF
            
        %Initial guesses for radii
        r1 = a_initial;
        r2 = a_final;
        
        %fmincon(@fSolveTofFunction, pi, [], [], [], [], safeTransferAngle, 2*pi + safeTransferAngle, [], opt_tf_angle);
        fminbnd(@fSolveTofFunction, safeTransferAngle, 2*pi + safeTransferAngle, opt_tf_angle);

        TOF_estimation = TOF_estimation * TOF_corrMult;
        tof_current = TOF_estimation;
    end

    if previousTime ~= currentTime
        T_nu = mod(currentTime - Tp1, P1);
        if T_nu ~= 0        
            nuHandle = @(angle) nuFromTime(angle, T_nu, n1, e1);

            nu1_i = fzero(nuHandle, [0, 2*pi], opt_nu_fzero);
        else
            nu1_i = 0;
        end
    
        theta1 = nu1_i + omega1;
        gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));
    
        r1 = p1 / (1+e1*cos(nu1_i));
        theta1_dot = sqrt(mju/a_initial^3) * a_initial^2/r1^2 * sqrt(1-e1^2);
        
        T_nu = mod(currentTime - Tp2, P2);
        if T_nu ~= 0
            nuHandle = @(angle) nuFromTime(angle, T_nu, n2, e2);
        
            nu2_i = fzero(nuHandle, [0, 2*pi], opt_nu_fzero);
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
    previousTime = currentTime;
    
    T_nu = mod(currentTime + tof_current - Tp2, P2);
    if T_nu ~= 0
        nuHandle = @(angle) nuFromTime(angle, T_nu, n2, e2);

        nu2 = fzero(nuHandle, [0, 2*pi], opt_nu_fzero);
    else
        nu2 = 0;
    end

    theta2 = nu2 + omega2;
    %theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
    
    theta_tilde = theta2 - theta1;
    if theta_tilde <= 0
        theta_tilde = theta_tilde + 2*pi;
    end

    gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
    r2 = p2 / (1+e2*cos(nu2));    
    r2_i = p2 / (1+e2*cos(nu2_i));    
    theta2_dot = sqrt(mju./a_final.^3).*a_final.^2./r2.^2.*sqrt(1-e2.^2);

    theta_f = 2.*pi.*N + theta_tilde;

    geometricAngle = acos(min(r1,r2)/max(r1,r2));

    if theta_f < geometricAngle
        theta_f = theta_f + 2*pi;
    end

    theta_vec = linspace(theta_0, theta_f, intApprox);
    theta_vec_acc = linspace(theta_0, theta_f, intApprox);

    paramVector = [mju, gamma1, gamma2, theta_f, theta1_dot, theta2_dot, r1, r2, theta1, theta2, nu2_i, r2_i];
    
%     gamma1
%     gamma2
%     theta_f
%     theta1_dot
%     theta2_dot
%     r1
%     r2
%     theta1
%     theta2
%     nu2_i
%     r2_i

    %% Check if TOF is a solution 

    unfeasibleOrbit = 0;

    %Check what the minimum d coefficient is
%     maxRadHandle = @(d_in) fMaxRadiusFunction(d_in, paramVector, rMax, a_initial);
% 
%     d_minimum = fzero(maxRadHandle, 0, opt_d_lim_fzero);

    %Limit the maximum distance geometrically
    if theta_f < pi

        objectDistance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta_f));
      
        alpha1 = acos((2*r1^2 - 2*r1*r2*cos(theta_f)) / (2*r1*objectDistance));

        objectToMax2 = sin(pi-alpha1) * objectDistance / sin(2*pi - theta_f);

        geometricMaxRadius = sqrt(objectToMax2^2 + r2^2);
        geometricMaxRadius = min(geometricMaxRadius, rMax);
    else
        geometricMaxRadius = rMax;
    end

    d_minimum = fFindRadiusFunction(paramVector, geometricMaxRadius);

%     tof_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0));
    %minFuncValue = min(fTimeMinReal(theta_vec_acc, d_minimum, paramVector, a_initial));
    
%     if ~isreal(tof_max)
        %Further refine the guess
%         opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
%         [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    
%         while ~isreal(tof_max)
        while min(fTimeMinReal(theta_vec_acc, d_minimum, paramVector, a_initial)) < 0
             %Move a relative amount towards larger values and an arbitary
             %small step to cross zero properly
            if d_minimum < 0
                d_minimum = d_minimum / dAdjustment + 1e-15;
            else
                d_minimum = d_minimum * dAdjustment + 1e-15;
            end
            %No real solutions exist for orbit
            if d_minimum > 1
                unfeasibleOrbit = 1;
                d_maximum = 1;
                return
            end
%             tof_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0));
            %minFuncValue = min(fTimeMinReal(theta_vec_acc, d_minimum, paramVector, a_initial));
%             [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
        end
%     end
%     %try
      %d_maximum = fzero(@fFindRadiusFunction, [d_minimum, 1], opt_d_lim_fzero);
      d_maximum = fFindRadiusFunction(paramVector, rMin);
%     %catch    
%         figure;
%         vectorD = linspace(d_minimum, d_minimum/1000, 1000);
%         plot(fFindRadiusFunction(vectorD));
%         
%         figure
%         plot(fInvRadiusZeroFunction(theta_f/2, vectorD));
%         
%         figure;
%     %end
    %d_maximum

    %tof_min = trapz(theta_vec_acc, fTimeMinReal(d_maximum, theta_vec_acc, paramVector, a_initial));

%     if ~isreal(tof_min)
%         opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
%         [~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle, d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    
        while ~isreal(fTimeMinReal(d_maximum, theta_vec_acc, paramVector, a_initial))
        %while minFuncValue < 0
            %Move a relative amount towards larger values and an arbitary
            %small step to cross zero properly
            if d_maximum < 0
                d_maximum = d_maximum * dAdjustment - 1e-15;
            else
                d_maximum = d_maximum / dAdjustment - 1e-15;
            end
            %tof_min = trapz(theta_vec_acc, fTimeMinReal(d_maximum, theta_vec_acc, paramVector, a_initial));
            %[~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle,
            %d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
            
%             x = cos(theta_vec_acc+theta1) .* fRadiusFunction(d_maximum, theta_vec_acc, 0);
%             y = sin(theta_vec_acc+theta1) .* fRadiusFunction(d_maximum, theta_vec_acc, 0);
%            
%             hold off;
%             plot(orbit1(1,:), orbit1(2,:), 'LineStyle',':', LineWidth=2);
%             hold on;
%             plot(orbit2(1,:), orbit2(2,:), 'LineStyle',':', LineWidth=2);
% 
%             plot(x, y, "Color", [0.2 0.7 0.2]);
            
        end
%     end
%     d_maximum
%     tof_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0))
%     tof_min2 = trapz(theta_vec_acc, fTimeFunction(d_maximum, theta_vec_acc, 0))

end

% %% Function for finding the largest radius of the transfer orbit
% % r = 1/f(theta), search roots of f(theta)
% function [inv_r_zero] = fInvRadiusZeroFunction(theta, d)
%     global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot
%     global a_initial; %For scaling the output
% 
%     a = 1/r1;
%     b = -tan(gamma1) / r1;
%     c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);
%     
%     efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
%                 -48*theta_f     18*theta_f^2 -2*theta_f^3; 
%                  20            -8*theta_f     theta_f^2];
%     
%     efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
%                 -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
%                 mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
%     
%     efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
%     
%     e = efg(1);
%     f = efg(2);
%     g = efg(3);
% 
%     inv_r_zero = a_initial * (a + b.*theta + c.*theta.^2 + d.*theta.^3 + e.*theta.^4 + f.*theta.^5 + g.*theta.^6);
% end

%% Function for finding the d-coefficient where orbit reaches minimum r
function [d_sol] = fFindRadiusFunction(paramVector, rTarget)
    %global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot
    %global rMin
    
    mju = paramVector(1);
    gamma1 = paramVector(2);
    gamma2 = paramVector(3);
    theta_f = paramVector(4);
    theta1_dot = paramVector(5);
    theta2_dot = paramVector(6);
    r1 = paramVector(7);
    r2 = paramVector(8);

%     a = 1/r1;
%     b = -tan(gamma1) / r1;
%     c = 1/(2*r1) * (mju / (r1^3 * theta1_dot^2) - 1);
%     
%     efg_Mat_1 = [30*theta_f^2  -10*theta_f^3  theta_f^4;
%                 -48*theta_f     18*theta_f^2 -2*theta_f^3; 
%                  20            -8*theta_f     theta_f^2];
%     
%     efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
%                 -tan(gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
%                 mju/(r2^4*theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
%     
%     efg = 1/(2*theta_f^6) * efg_Mat_1 * efg_Mat_2;
%     
%     e = efg(1);
%     f = efg(2);
%     g = efg(3);
% 
%     r_error = 1./(a + b.*(theta_f/2) + c.*(theta_f/2).^2 + d.*(theta_f/2).^3 + e.*(theta_f/2).^4 + f.*(theta_f/2).^5 + g.*(theta_f/2).^6) - rMin;


d_sol = (64*(1/rTarget - 1/r1 + (theta_f^6*(((5*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (10*tan(gamma1)*theta_f)/r1 + 10/r1 - 10/r2)/theta_f^6 - ((4*tan(gamma2))/r2 - (4*tan(gamma1))/r1 + (4*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^5 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^4)))/64 + (theta_f^4*(((15*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/(2*r1) - (15*tan(gamma1)*theta_f)/r1 + 15/r1 - 15/r2)/theta_f^4 - ((5*tan(gamma2))/r2 - (5*tan(gamma1))/r1 + (5*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^3 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/(2*theta_f^2)))/16 - (theta_f^5*(((12*(mju/(r1^3*theta1_dot^2) - 1)*theta_f^2)/r1 - (24*tan(gamma1)*theta_f)/r1 + 24/r1 - 24/r2)/theta_f^5 - ((9*tan(gamma2))/r2 - (9*tan(gamma1))/r1 + (9*theta_f*(mju/(r1^3*theta1_dot^2) - 1))/r1)/theta_f^4 + ((mju/(r1^3*theta1_dot^2) - 1)/r1 + 1/r2 - mju/(r2^4*theta2_dot^2))/theta_f^3))/32 + (theta_f*tan(gamma1))/(2*r1) - (theta_f^2*(mju/(r1^3*theta1_dot^2) - 1))/(8*r1)))/theta_f^3;
 
end
%% Function for finding the initial d-coefficient form furthest allowed point 
function [r_error] = fMaxRadiusFunction(d, paramVector, rMax, a_initial)
       
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
    
    r_error = a_initial.*((a + b.*(theta_f/2) + c.*(theta_f/2).^2 + d.*(theta_f/2).^3 + e.*(theta_f/2).^4 + f.*(theta_f/2).^5 + g.*(theta_f/2).^6) - 1./rMax);
     
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

    %global a_initial %For scaling

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

    %     timeTroublePart = (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4);
    % 
    %     timeTroubleDiff = diff(timeTroublePart);
    % 
    %     timeTroubleZero = solve(timeTroubleDiff == 0, theta);
    % 
    %     timeTroubleRoot = vpa(timeTroubleZero)
    % 


end

%% For Solving the initial TOF guess 
function [meetAngleError] = fSolveTofFunction(transferAngle)
    global currentTime N mju TOF_estimation
    global n1 e1 omega1 P1 p1 r1 Tp1
    global n2 e2 omega2 P2 p2 r2 Tp2
    %global TOF_corrMult;
    global n_nu e_nu T_nu
    global theta1 theta2 theta_f

    global opt_nu_fzero

    TOF_estimation = (2*N*pi + transferAngle)*sqrt((r1+r2)^3/(8*mju));
    %TOF_estimation = (2*N*pi + modifiedAngle)*sqrt((r1+r2)^3/(8*mju));
   
    T_nu = mod(currentTime - Tp1, P1);
    if T_nu ~= 0
        %Used in the next fzero function
        n_nu = n1;
        e_nu = e1;

        nuHandle = @(angle) nuFromTime(angle, T_nu, n_nu, e_nu);

        nu1_i = fzero(nuHandle, [0, 2*pi], opt_nu_fzero);
    else
        nu1_i = 0;
    end

    theta1 = nu1_i + omega1;

    r1 = p1 / (1+e1*cos(nu1_i));

    %Used in the next fzero function
    T_nu = mod(currentTime + TOF_estimation - Tp2, P2);
    if T_nu ~= 0

        n_nu = n2;
        e_nu = e2;    
    
        nuHandle = @(angle) nuFromTime(angle, T_nu, n_nu, e_nu);

        nu2 = fzero(nuHandle, [0, 2*pi], opt_nu_fzero);
    else
        nu2 = 0;
    end

    theta2 = nu2 + omega2;
    theta_tilde = theta2 - theta1;
        if theta_tilde < 0
            theta_tilde = theta_tilde + 2*pi;
        end
    %theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
%     theta_tilde_optim = theta_tilde
    
    r2 = p2 / (1+e2*cos(nu2));    
    
    theta_f = 2.*pi.*N + theta_tilde;

%     transferAngle;
%     theta_f2 = theta_f;
%     TOF_estimation;
    meetAngleError = transferAngle - theta_f;
end

% function [trueAnomaly] = fNuApproxFunction(meanAnomaly, eccentricity)
% 
%     if eccentricity == 0
%         trueAnomaly = meanAnomaly;
%         return
%     end
%     
%     if meanAnomaly > pi
%         ratio = 1-(meanAnomaly-pi)./pi;
%     else
%         ratio = meanAnomaly./pi;
%     end
% 
%     c = 2-1/sqrt(eccentricity);
% 
%     a = 1/2 * (sqrt(5-4*c)-1);
% 
%     trueAnomaly = pi.*(1-((1-c)./(ratio+a)-a));
% 
%     if meanAnomaly > pi
%         trueAnomaly = 2*pi - trueAnomaly;
%     end
% 
% end