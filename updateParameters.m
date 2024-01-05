function [] = updateParameters(updateTOF)
    % Calculating the orbital parameters
    %initial orbit parameters
    global a_initial a_final currentTime mju N tof_current;
    global theta1 theta2 theta_f omega1 omega2 e1 e2 theta2_dot;
    global gamma1 r1 gamma2 r2 theta1_dot theta_vec;
    global nu1_i nu2_i r1_i r2_i intApprox;
    global Tp1 Tp2 TOF_estimation;
    global d_minimum d_maximum;% rMin rMax;
    global theta_0 e_nu n_nu T_nu;
    global n1 P1 n2 P2 p1 p2;
    global previousTime unfeasibleOrbit safeTransferAngle;
    global TOF_corrMult dAdjustment;

    global theta_vec_acc;
    
    global opt_nu_fzero opt_tf_angle opt_d_lim_fzero;

    if updateTOF
            
        %Initial guesses for radii
        r1 = a_initial;
        r2 = a_final;
        
        %fmincon(@fSolveTofFunction, pi, [], [], [], [], safeTransferAngle, 2*pi + safeTransferAngle, [], opt_tf_angle);
        fminbnd(@fSolveTofFunction, safeTransferAngle, 2*pi + safeTransferAngle, opt_tf_angle);

%         bestAngle = 0;
%         minimumError = Inf;
% 
%         angleVec = linspace(safeTransferAngle, (1+N)*2*pi + safeTransferAngle, 256);
%         valVec = zeros(256,1);
% 
%         for i = 1:256
%             valVec(i) = fSolveTofFunction(angleVec(i));
%         end
% 
%         figure;
%         plot(angleVec, valVec)
% 
%         figure;
%         for tfAngleGuess = linspace(safeTransferAngle, (1+N)*2*pi + safeTransferAngle, 24)
%             angleError = abs(fSolveTofFunction(tfAngleGuess));
%             if angleError < minimumError
%                 minimumError  = angleError;
%                 bestAngle = tfAngleGuess;
%             end
%         end
%         
%         fSolveTofFunction(bestAngle);

%         %Mean anomaly
%         syms tof
% 
%         m1 = n1 * mod(currentTime - Tp1, P1);
%         m2 = n1 * mod(currentTime + tof - Tp2, P2);
% 
%         if m1 < pi
%             nuApprox1 = fNuApproxFunction(m1, e1)
% 
% 
%         tfAngle = mod(m2-m1 + 2*pi, 2*pi);

%         time_vec = linspace(0, P1-1, 1000);
% 
%         mean_vec = n1 * mod(time_vec - 0, P1);
% 
%         figure;
%         hold on;
%         correctVec = zeros(256,1);
%         for e = linspace(0, 0.999, 20)
%             for i = 1:1000
%                 correctVec(i) = fNuApproxFunction(mean_vec(i), e);
%             end
%             plot(mean_vec, correctVec);
%         end

        %tfAngleApprox = pi;
% 
%         r1_approx = a_initial*(1-0.5*e1);
%         r2_approx = a_final*(1-0.5*e2);

%         tof_guess = (pi + 2*pi*N)*sqrt((r1_approx+r2_approx)^3/(8*mju));
% 
%         m1 = n1 * mod(currentTime - Tp1, P1);
%         m2 = n1 * mod(currentTime + tof_guess - Tp2, P2);
%         
%         nuApprox1 = fNuApproxFunction(m1, e1);
%         nuApprox2 = fNuApproxFunction(m2, e2);
% 
%         tfAngleApprox = mod(nuApprox2 + omega2 - (nuApprox1 + omega1) + 2*pi, 2*pi);
% 
%         r1_approx = p1/ (1+e1*cos(nuApprox1));
%         r2_approx = p2/ (1+e2*cos(nuApprox2));
%         
%         tof_approx = (tfAngleApprox + 2*pi*N)*sqrt((r1_approx+r2_approx)^3/(8*mju));
% 
%         figure;
%         hold on;
% 
%         bestTOFApprox = 0;
%         minTOFError = Inf;
% 
%         bestTfAngle = 0;
%         tof_approx = 0;
% 
%         for tfAngleGuess = linspace(pi/2, 2*pi+pi/2, 1000)
%     
%             m1 = n1 * mod(currentTime - Tp1, P1);
%             m2 = n1 * mod(currentTime + tof_approx - Tp2, P2);
%             
%             nuApprox1 = fNuApproxFunction(m1, e1);
%             nuApprox2 = fNuApproxFunction(m2, e2);
%     
%             tfAngleApprox = mod(nuApprox2 + omega2 - (nuApprox1 + omega1) + 2*pi, 2*pi);
%             
%             if tfAngleApprox < safeTransferAngle
%                 tfAngleApprox = tfAngleApprox + 2*pi;
%             end
% 
%             r1_approx = p1/ (1+e1*cos(nuApprox1));
%             r2_approx = p2/ (1+e2*cos(nuApprox2));
%             
%             tof_approx_new = (tfAngleApprox + 2*pi*N)*sqrt((r1_approx+r2_approx)^3/(8*mju));
% 
%             tofError = abs(tof_approx_new - tof_approx);
% 
%             tof_approx = tof_approx_new;
% 
%             if tofError < minTOFError
%                 bestTfAngle = tfAngleApprox;
%                 minTOFError = tofError;
%                 bestTOFApprox = tof_approx;
%             end
% 
%             plot(tfAngleGuess, abs(tfAngleApprox - tfAngleGuess) ,'or', 'MarkerSize',5,'MarkerFaceColor','g');
%             %plot(tfAngleGuess, r1_approx,'or', 'MarkerSize',5,'MarkerFaceColor','g');
%             %plot(tfAngleGuess, r2_approx,'or', 'MarkerSize',5,'MarkerFaceColor','b');
%             %plot(tfAngleGuess, tof_approx ,'or', 'MarkerSize',5,'MarkerFaceColor','g')
%             %plot(tfAngleGuess, bestTOFApprox,'or', 'MarkerSize',5,'MarkerFaceColor','r')
%         end
% 
         %TOF_estimation1 = bestTOFApprox;%TOF_corrMult*(1+2*N)*pi*sqrt((a_initial*(1-0.5*e1)+a_final*(1-0.5*e2))^3/(8*mju));
         %TOF_estimation1 = TOF_corrMult*(1+2*N)*pi*sqrt((a_initial*(1-0.5*e1)+a_final*(1-0.5*e2))^3/(8*mju));
         %TOF_estimation2 = TOF_corrMult*(1+2*N)*pi*sqrt((a_initial+a_final)^3/(8*mju));

         TOF_estimation = TOF_estimation * TOF_corrMult;

         tof_current = TOF_estimation;
    end

    if previousTime ~= currentTime
        T_nu = mod(currentTime - Tp1, P1);
        if T_nu ~= 0
            %Used in the next fzero function
            n_nu = n1;
            e_nu = e1;
        
            nu1_i = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
        else
            nu1_i = 0;
        end
    
        theta1 = nu1_i + omega1;
        gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));
    
        r1 = p1 / (1+e1*cos(nu1_i));
        r1_i = p1 / (1+e1*cos(nu1_i));
        theta1_dot = sqrt(mju/a_initial^3) * a_initial^2/r1^2 * sqrt(1-e1^2);
        
        T_nu = mod(currentTime - Tp2, P2);
        if T_nu ~= 0
            %Used in the next fzero function
            n_nu = n2;
            e_nu = e2;
        
            nu2_i = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
        else
            nu2_i = 0;
        end
    end
    previousTime = currentTime;
    
    T_nu = mod(currentTime + tof_current - Tp2, P2);
    if T_nu ~= 0
        %Used in the next fzero function
        n_nu = n2;
        e_nu = e2;
        
        nu2 = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
    else
        nu2 = 0;
    end

    theta2 = nu2 + omega2;
    theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);

%     theta_tilde_resolve = theta_tilde
    
    gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
    r2 = p2 / (1+e2*cos(nu2));    
    r2_i = p2 / (1+e2*cos(nu2_i));    
    theta2_dot = sqrt(mju./a_final.^3).*a_final.^2./r2.^2.*sqrt(1-e2.^2);

%     theta1
%     theta2
    theta_f = 2.*pi.*N + theta_tilde;

    if theta_f < safeTransferAngle
        theta_f = theta_f + 2*pi;
    end

    theta_vec = linspace(theta_0, theta_f, intApprox);
    theta_vec_acc = linspace(theta_0, theta_f, 1000);

    %% Check if TOF is a solution 

    unfeasibleOrbit = 0;

    %Check what the minimum d coefficient is

%     answer0 = fMaxRadiusFunction(0)
%     answer1 = fMaxRadiusFunction(-1e-12)
%     answer2 = fMaxRadiusFunction(1e-12)
%     

    d_minimum = fzero(@fMaxRadiusFunction, 0, opt_d_lim_fzero);

%     tof_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0));
    minFuncValue = min(fTimeMinReal(theta_vec, d_minimum));
    
%     if ~isreal(tof_max)
        %Further refine the guess
%         opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
%         [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    
%         while ~isreal(tof_max)
        while minFuncValue < 0
             %Move a relative amount towards larger values and an arbitary
             %small step to cross zero properly
            if d_minimum < 0
                d_minimum = d_minimum / dAdjustment + 1e-18;
            else
                d_minimum = d_minimum * dAdjustment + 1e-18;
            end
            %No real solutions exist for orbit
            if d_minimum > 1
                unfeasibleOrbit = 1;
                d_maximum = 1;
                return
            end
%             tof_max = trapz(theta_vec, fTimeFunction(d_minimum, theta_vec, 0));
            minFuncValue = min(fTimeMinReal(theta_vec, d_minimum));
%             [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
        end
%     end
%     %try
      d_maximum = fzero(@fMinRadiusFunction, [d_minimum, 1], opt_d_lim_fzero);
%     %catch    
%         figure;
%         vectorD = linspace(d_minimum, d_minimum/1000, 1000);
%         plot(fMinRadiusFunction(vectorD));
%         
%         figure
%         plot(fInvRadiusZeroFunction(theta_f/2, vectorD));
%         
%         figure;
%     %end
    %d_maximum

    tof_min = trapz(theta_vec_acc, fTimeFunction(d_maximum, theta_vec_acc, 0));

%     if ~isreal(tof_min)
%         opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
%         [~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle, d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    
        while ~isreal(tof_min)
        %while minFuncValue < 0
            %Move a relative amount towards larger values and an arbitary
            %small step to cross zero properly
            if d_maximum < 0
                d_maximum = d_maximum * dAdjustment - 1e-18;
            else
                d_maximum = d_maximum / dAdjustment - 1e-18;
            end
            tof_min = trapz(theta_vec_acc, fTimeFunction(d_maximum, theta_vec_acc, 0));
            %[~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle, d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
        end
%     end
%     d_maximum
%     tof_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0))
%     tof_min2 = trapz(theta_vec_acc, fTimeFunction(d_maximum, theta_vec_acc, 0))

end

%% Function for finding the largest radius of the transfer orbit
% r = 1/f(theta), search roots of f(theta)
function [inv_r_zero] = fInvRadiusZeroFunction(theta, d)
    global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot
    global a_initial; %For scaling the output

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

    inv_r_zero = a_initial * (a + b.*theta + c.*theta.^2 + d.*theta.^3 + e.*theta.^4 + f.*theta.^5 + g.*theta.^6);
end

%% Function for finding the d-coefficient where orbit reaches minimum r
function [r_error] = fMinRadiusFunction(d)
    global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot
    global rMin

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

    r_error = 1./(a + b.*(theta_f/2) + c.*(theta_f/2).^2 + d.*(theta_f/2).^3 + e.*(theta_f/2).^4 + f.*(theta_f/2).^5 + g.*(theta_f/2).^6) - rMin;
end
%% Function for finding the initial d-coefficient form furthest allowed point 
function [r_error] = fMaxRadiusFunction(d)
    
    global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot;
    global rMax a_initial;
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
function [angleError] = nuFromTime(nu_time)
    global e_nu n_nu T_nu
    %nu_time
    if nu_time <= pi
        E = 2*atan(tan(nu_time/2)/sqrt((1+e_nu)/(1-e_nu)));
    else
        E = 2*pi + 2*atan(tan(nu_time/2)/sqrt((1+e_nu)/(1-e_nu)));
    end

    angleError = E-e_nu*sin(E) - n_nu*T_nu;
end

%% Solve maximum d-coefficient from part of the time function
function [timeImgPart] = fTimeMinReal(theta, d)
    global mju gamma1 gamma2 theta_f theta1_dot theta2_dot r1 r2;
    global a_initial %For scaling

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
function [meetAngleError] = fSolveTofFunction(transferAngle)
    global currentTime N mju TOF_estimation
    global n1 e1 omega1 P1 p1 r1 Tp1
    global n2 e2 omega2 P2 p2 r2 Tp2
    global TOF_corrMult;
    global n_nu e_nu T_nu
    global theta1 theta2 theta_f

    global opt_nu_fzero
    
    transferAngle;
    TOF_estimation = (2*N*pi + transferAngle)*sqrt((r1+r2)^3/(8*mju));
    
    T_nu = mod(currentTime - Tp1, P1);
    if T_nu ~= 0
        %Used in the next fzero function
        n_nu = n1;
        e_nu = e1;
    
        nu1_i = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
    else
        nu1_i = 0;
    end

    theta1 = nu1_i + omega1;

    r1 = p1 / (1+e1*cos(nu1_i));

    %Used in the next fzero function
    n_nu = n2;
    e_nu = e2;
    T_nu = mod(currentTime + TOF_estimation - Tp2, P2);
    
    nu2 = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);

    theta2 = nu2 + omega2;
    theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
%     theta_tilde_optim = theta_tilde
    
    r2 = p2 / (1+e2*cos(nu2));    
    
    theta_f = 2.*pi.*N + theta_tilde;

%     transferAngle;
%     theta_f2 = theta_f;
%     TOF_estimation;
    meetAngleError = transferAngle - theta_f;
end

function [trueAnomaly] = fNuApproxFunction(meanAnomaly, eccentricity)

    if eccentricity == 0
        trueAnomaly = meanAnomaly;
        return
    end
    
    if meanAnomaly > pi
        ratio = 1-(meanAnomaly-pi)./pi;
    else
        ratio = meanAnomaly./pi;
    end

    c = 2-1/sqrt(eccentricity);

    a = 1/2 * (sqrt(5-4*c)-1);

    trueAnomaly = pi.*(1-((1-c)./(ratio+a)-a));

    if meanAnomaly > pi
        trueAnomaly = 2*pi - trueAnomaly;
    end

end