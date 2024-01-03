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
    global previousTime unfeasibleOrbit safeTransferAngle TOF_corrMult;
    

    if updateTOF
            
%         %Initial guesses for radii
%         r1 = a_initial;
%         r2 = a_final;
%         
% 
%         %opt_tf_angle = optimset('TolFun', 1e-3);
%         %fzero(@fSolveTofFunction, [safeTransferAngle, (1+N)*2*pi + safeTransferAngle], opt_tf_angle);
% 
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

        TOF_estimation = TOF_corrMult*(1+2*N)*pi*sqrt((a_initial*(1-0.5*e1)+a_final*(1-0.5*e2))^3/(8*mju));
        
        tof_current = TOF_estimation;
    end

    opt_nu_fzero = optimset('TolFun', 1e-3);

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

    %Used in the next fzero function
    n_nu = n2;
    e_nu = e2;
    T_nu = mod(currentTime + tof_current - Tp2, P2);
    
    nu2 = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);

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

    %% Check if TOF is a solution 

    unfeasibleOrbit = 0;

    %Check what the minimum d coefficient is
    opt_fzero = optimset('TolFun', 1e1, 'TolX', 1e-15, 'Display', 'off');

%     answer0 = fMaxRadiusFunction(0)
%     answer1 = fMaxRadiusFunction(-1e-12)
%     answer2 = fMaxRadiusFunction(1e-12)
%     

    d_minimum = fzero(@fMaxRadiusFunction, 0, opt_fzero);

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
                d_minimum = d_minimum / 1.1 + 1e-15;
            else
                d_minimum = d_minimum * 1.1 + 1e-15;
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

    d_maximum = fzero(@fMinRadiusFunction, [d_minimum, 1], opt_fzero);

    %d_maximum

    tof_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0));

%     if ~isreal(tof_min)
%         opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-4, 'Display', 'off');
%         [~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle, d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    
        while ~isreal(tof_min)
        %while minFuncValue < 0
            %Move a relative amount towards larger values and an arbitary
            %small step to cross zero properly
            if d_maximum < 0
                d_maximum = d_maximum * 1.01 - 1e-15;
            else
                d_maximum = d_maximum / 1.01 - 1e-15;
            end
            tof_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0));
            %[~, minFuncValue] = fmincon(@(angle) fTimeMinReal(angle, d_maximum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
        end
%     end
    %d_maximum
    %tof_min = trapz(theta_vec, fTimeFunction(d_maximum, theta_vec, 0))

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
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e_nu)/(1-e_nu)));
    angleError = pi+E-e_nu*sin(E) - n_nu*T_nu;
    %angleError = mod(angleError)
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
    
    transferAngle;
    TOF_estimation = (2*N*pi + transferAngle)*sqrt((r1+r2)^3/(8*mju));

    opt_nu_fzero = optimset('TolFun', 1e-3);
    
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