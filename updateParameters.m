function [] = updateParameters(updateTOF)
    % Calculating the orbital parameters
    %initial orbit parameters
    global a_initial a_final currentTime mju N tof_current;
    global theta1 theta2 theta_f omega1 omega2 e1 e2 theta2_dot;
    global gamma1 r1 P1 gamma2 r2 P2 theta1_dot theta_vec;
    global nu1_i nu2_i r1_i r2_i intApprox;
    global Tp1 Tp2 TOF_estimation;
    global d_minimum d_maximum rMin rMax;
    global theta_0 timeFunction_nn;

    if updateTOF
        TOF_estimation = (1+2*N)*pi*sqrt((a_initial+a_final)^3/(8*mju));
        tof_current = TOF_estimation;
    end

    a1 = a_initial;
    n1 = sqrt(mju/a1^3);
    P1 = 2*pi/n1;
    T1_i = mod(currentTime - Tp1, P1);

    syms nu_time
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e1)/(1-e1)));
    nuSolver = n1*T1_i==pi+E-e1*sin(E);
    nuSolutions_i = vpasolve(nuSolver, nu_time);
    nuSolutions_i = double(nuSolutions_i);
    nu1_i = mod(nuSolutions_i + 2*pi, 2*pi);
    
    theta1 = nu1_i + omega1;
    gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));

    p1 = a1 * (1-e1^2);
    r1 = p1 / (1+e1*cos(nu1_i));
    r1_i = p1 / (1+e1*cos(nu1_i));
    theta1_dot = sqrt(mju/a1^3) * a1^2/r1^2 * sqrt(1-e1^2);
    
    a2 = a_final;
    n2 = sqrt(mju/a2^3);
    P2 = 2*pi/n2;
    T2_i = mod(currentTime - Tp2, P2);
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
    nuSolver = n2*T2_i==pi+E-e2*sin(E);
    nuSolutions_i = vpasolve(nuSolver, nu_time);
    nuSolutions_i = double(nuSolutions_i);
    nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);

    finalTime = currentTime + tof_current;
    T = mod(finalTime - Tp2, P2);
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
    nuSolver = n2*T==pi+E-e2*sin(E);
    nuSolutions_f = vpasolve(nuSolver, nu_time);
    nuSolutions_f = double(nuSolutions_f);
    
    nu2 = mod(nuSolutions_f + 2*pi, 2*pi);
    theta2 = nu2 + omega2;
    %theta_tilde = theta2 - theta1;
    theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
    
    gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
    p2 = a2 * (1-e2^2);
    r2 = p2 / (1+e2*cos(nu2));    
    r2_i = p2 / (1+e2*cos(nu2_i));    
    theta2_dot = sqrt(mju./a2.^3).*a2.^2./r2.^2.*sqrt(1-e2.^2);
    
    theta_f = 2.*pi.*N + theta_tilde;

%     safeTransferAngle = pi/2;
%     if theta_f < safeTransferAngle
%         theta_f = theta_f + 2*pi;
%     end

    %% Solving the coefficients
    syms d theta;

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
    gamma = atan(-r * (b + 2*c*theta + 3*d*theta^2 + 4*e*theta^3 + 5*f*theta^4 + 6*g*theta^5));
    
    global thrustFunction thetaDotFunction thetaDotSquareFunction timeFunction radiusFunction;

    thetaDotFunction = sqrt((mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
    thrustFunction = -mju / (2 * r^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta^2 + 120*g*theta^3 - tan(gamma)/r) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)^2;
    
    timeFunction = sqrt((r^4/mju) * (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));

    timeFunction_nn = @(d_coeff, angle) double(subs(subs(timeFunction, d, d_coeff), theta, angle));

    thetaDotSquareFunction = (mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4);

    radiusFunction = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);

    theta_vec = linspace(theta_0, theta_f, intApprox);


    %% Check if TOF is a solution 
% 
%     %Check what the minimum d coefficient is
%     radiusMax_n = subs(1/radiusFunction, theta, theta_f/2);
%     radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
%     
%     %Move the d-coefficient to the positive direction -->
%     %Reduce chance to find orbits going into infinity
%     safe_d_multiplier = 2;%1.4;
%     d_minimum = fzero(radiusMax_nn, 0);

    radiusMax_n = subs(1e9/radiusFunction - 1e9/rMax, theta, theta_f/2);
    radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
     
    opt = optimset('TolFun', 1e-15, 'TolX', 1e-15, 'Display', 'iter');
    d_minimum = fzero(radiusMax_nn, 0, opt);

    invR_Function = @(angle) double(subs(subs(1e9/radiusFunction, d, d_minimum), theta, angle));
    [minPoints, minFuncValue,EXITFLAG,OUTPUT,LAMBDA] = fmincon(invR_Function, theta_f/2, [], [], [], [], theta_0, theta_f, [], opt);
    
%     minFuncValue
%     minPoints
%     d_minimum
%     EXITFLAG
%     OUTPUT.message
%     LAMBDA
%     abaddu=1/radiusFunction
    
%     figure;
%     plot(theta_vec, invR_Function(theta_vec))

%ADD CHECK FOR ZERO BEING A WRONG ANSWER


    while minFuncValue < 0
        if d_minimum < 0
            d_minimum = d_minimum / 2;
        else
            d_minimum = d_minimum * 2;
        end
        invR_Function = @(angle) double(subs(subs(1e9/radiusFunction, d, d_minimum), theta, angle));
        [minPoints, minFuncValue,EXITFLAG,OUTPUT,LAMBDA] = fmincon(invR_Function, theta_f/2, [], [], [], [], theta_0, theta_f, [], opt);
%         minFuncValue
%         minPoints
%         d_minimum
%         EXITFLAG
%         OUTPUT.message
%         LAMBDA
    end


%     if d_minimum < 0
%         d_minimum = d_minimum / safe_d_multiplier;
%     else
%         d_minimum = d_minimum * safe_d_multiplier;
%     end
    %d_minimum = d_minimum + abs(d_minimum) * safe_d_multiplier;
    

    %Could additionally be defined as a maximum distance as with the upper
    %limit of the d-coefficient
   
%     d_minimum
%     radiusFunction
%     theta_f
%     rMin

    radiusMin_n = subs(radiusFunction - rMin, theta, theta_f/2);
    radiusMin_nn = @(d_max_in) double(subs(radiusMin_n, d, d_max_in));
    
    d_maximum = fzero(radiusMin_nn, [d_minimum, 1]);

end
