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
    global previousTime;

    if updateTOF
        TOF_estimation = (1+2*N)*pi*sqrt((a_initial*(1-0.5*e1)+a_final*(1-0.5*e2))^3/(8*mju));
        tof_current = TOF_estimation;
    end

    opt_nu_fzero = optimset('TolFun', 1e-3);
    
    if previousTime ~= currentTime
        if Tp1 ~= 0
            %Used in the next fzero function
            n_nu = n1;
            e_nu = e1;
            T_nu = mod(currentTime - Tp1, P1);
        
            nu1_i = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
            %nu1_i = mod(nu1_i + 2*pi, 2*pi);
        else
            nu1_i = 0;
        end
    
        theta1 = nu1_i + omega1;
        gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));
    
        r1 = p1 / (1+e1*cos(nu1_i));
        r1_i = p1 / (1+e1*cos(nu1_i));
        theta1_dot = sqrt(mju/a_initial^3) * a_initial^2/r1^2 * sqrt(1-e1^2);
    
        if Tp2 ~= 0
            %Used in the next fzero function
            n_nu = n2;
            e_nu = e2;
            T_nu = mod(currentTime - Tp2, P2);
        
            nu2_i = fzero(@nuFromTime, [0, 2*pi], opt_nu_fzero);
            %nu2_i = mod(nu2_i + 2*pi, 2*pi);
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
    %nu2 = mod(nu2 + 2*pi, 2*pi);

    theta2 = nu2 + omega2;
    theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
    
    gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
    r2 = p2 / (1+e2*cos(nu2));    
    r2_i = p2 / (1+e2*cos(nu2_i));    
    theta2_dot = sqrt(mju./a_final.^3).*a_final.^2./r2.^2.*sqrt(1-e2.^2);
    
    theta_f = 2.*pi.*N + theta_tilde;

    safeTransferAngle = pi;
    if theta_f < safeTransferAngle
        theta_f = theta_f + 2*pi;
    end


    %% Solving the coefficients
    syms d theta;
% 
%     mju
%     r1
%     r2
%     gamma1
%     gamma2
%     theta_f
%     theta1_dot
%     theta2_dot

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
    
    r = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);
    
    global thetaDotFunction thetaDotSquareFunction radiusFunction;

    thetaDotFunction = sqrt((mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
    
    thetaDotSquareFunction = (mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4);

    radiusFunction = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);

    theta_vec = linspace(theta_0, theta_f, intApprox);

    %% Check if TOF is a solution 

%     %Check what the minimum d coefficient is
%     radiusMax_n = subs(1/radiusFunction, theta, theta_f/2);
%     radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
%     
%     %Move the d-coefficient to the positive direction -->
%     %Reduce chance to find orbits going into infinity
%     safe_d_multiplier = 2;%1.4;
%     d_minimum = fzero(radiusMax_nn, 0);

%     radiusMax_n = subs(1e9/radiusFunction - 1e9/rMax, theta, theta_f/2);
%     radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
     
    opt_fzero = optimset('TolFun', 1e1, 'TolX', 1e-15, 'Display', 'off');
    d_minimum = fzero(@fMaxRadiusFunction, 0, opt_fzero);
    
    opt_fmincon = optimset('TolFun', 1e-3, 'TolX', 1e-15, 'Display', 'off');
    [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);


    while minFuncValue < 0
         %Move a relative amount towards larger values and an arbitary
         %small step to cross zero properly
        if d_minimum < 0
            d_minimum = d_minimum / 2 + 1e-18;
        else
            d_minimum = d_minimum * 2 + 1e-18;
        end

        [~, minFuncValue] = fmincon(@(angle) fInvRadiusZeroFunction(angle, d_minimum), theta_f/2, [], [], [], [], theta_0, theta_f, [], opt_fmincon);
    end
 
    d_maximum = fzero(@fMinRadiusFunction, [d_minimum, 1], opt_fzero);
    
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

