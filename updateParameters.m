function [] = updateParameters(updateTOF)
    % Calculating the orbital parameters
    %initial orbit parameters
    global a_initial a_final currentTime mju N tof_current;
    global theta1 theta2 theta_f omega1 omega2 e1 e2 theta2_dot;
    global gamma1 r1 P1 gamma2 r2 P2 theta1_dot;
    global nu1_i nu2_i r1_i r2_i;
    global Tp1 Tp2 TOF_estimation;
    global d_minimum d_maximum rMin;

    if updateTOF
        TOF_estimation = (1+2*N)*pi*sqrt((a_initial+a_final)^3/(8*mju));

        tof_current = TOF_estimation;
    end

    %tf = 86400 * 365.25 * 0.5;
    %ratio = tf/TOF_estimation;

    a1 = a_initial;
    %calulating new point for object 1 as well
    %Time of last perigee pass for object 1
    %Tp1 = 0;
    
    %Trying this one instead
    n1 = sqrt(mju/a1^3);
    P1 = 2*pi/n1;
    T1_i = mod(currentTime - Tp1, P1);
    %T1_i = currentTime - Tp1;
    
    % nT = E - e*sin(E)
    % E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));
    
    syms nu_time
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e1)/(1-e1)));
    nuSolver = n1*T1_i==pi+E-e1*sin(E);
    nuSolutions_i = vpasolve(nuSolver, nu_time);
    nuSolutions_i = double(nuSolutions_i);
    
    nu1_i = mod(nuSolutions_i + 2*pi, 2*pi);
    
    %Now we can continue by calculating the true anomaly when the spacecraft
    %reaches orbit 2. At currentTime + tf
    
%     T1 = mod(currentTime, P1) - Tp1;
%     E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e1)/(1-e1)));
%     nuSolver = n1*T1==pi+E-e1*sin(E);
%     nuSolutions_f = vpasolve(nuSolver, nu_time);
%     nuSolutions_f = double(nuSolutions_f);
%     
%     nu1 = mod(nuSolutions_f + 2*pi, 2*pi);

    %theta1 = nu1 - omega1;
    theta1 = nu1_i - omega1;
    %gamma1 = asin(e1 * sin(nu1) / sqrt(1+2*e1*cos(nu1) + e1^2));
    gamma1 = asin(e1 * sin(nu1_i) / sqrt(1+2*e1*cos(nu1_i) + e1^2));
    p1 = a1 * (1-e1^2);
    %r1 = p1 / (1+e1*cos(nu1));
    r1 = p1 / (1+e1*cos(nu1_i));
    r1_i = p1 / (1+e1*cos(nu1_i));
    theta1_dot = sqrt(mju/a1^3) * a1^2/r1^2 * sqrt(1-e1^2);
    

    %Target orbit parameters
    %Some are known, while others are calculated from desired tof and initial
    %conditions
    a2 = a_final;
    %Time of last perigee pass for object 2
    %Tp2 = 0;
    
    %Trying this one instead
    n2 = sqrt(mju/a2^3);
    P2 = 2*pi/n2;
    T2_i = mod(currentTime - Tp2, P2);
    %T2_i = currentTime - Tp2;

    syms nu_time
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
    nuSolver = n2*T2_i==pi+E-e2*sin(E);
    nuSolutions_i = vpasolve(nuSolver, nu_time);
    nuSolutions_i = double(nuSolutions_i);
    
    nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);

    %Now we can continue by calculating the true anomaly when the spacecraft
    %reaches orbit 2. At currentTime + tof_current
    
    finalTime = currentTime + tof_current;
    T = mod(finalTime - Tp2, P2);
    %T = finalTime - Tp2;

    % nT = E - e*sin(E)
    % E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));
    
    syms nu_time
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
    nuSolver = n2*T==pi+E-e2*sin(E);
    nuSolutions_f = vpasolve(nuSolver, nu_time);
    nuSolutions_f = double(nuSolutions_f);
    
    nu2 = mod(nuSolutions_f + 2*pi, 2*pi);
    %In reference coords
    theta2 = nu2 - omega2;
    theta_tilde = mod(theta2 - theta1 + 2*pi, 2*pi);
    %theta_tilde = theta2 - theta1;
    
    gamma2 = asin(e2 * sin(nu2) / sqrt(1+2*e2*cos(nu2) + e2^2));
    p2 = a2 * (1-e2^2);
    r2 = p2 / (1+e2*cos(nu2));    
    r2_i = p2 / (1+e2*cos(nu2_i));    
    theta2_dot = sqrt(mju./a2.^3).*a2.^2./r2.^2.*sqrt(1-e2.^2);
    
    %The total transfer angle is represented by:
    
    theta_f = 2.*pi.*N + theta_tilde;

    % Solving the coefficients
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

    thetaDotSquareFunction = (mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4);

    radiusFunction = 1 / (a + b*theta + c*theta^2 + d*theta^3 + e*theta^4 + f*theta^5 + g*theta^6);

    %% Check if TOF is a solution 

    %Check what the minimum d coefficient is
    radiusMax_n = subs(1/radiusFunction, theta, theta_f/2);
    radiusMax_nn = @(d_min_in) double(subs(radiusMax_n, d, d_min_in));
    
    %Needs to be <1
    safe_d_multiplier = 0.9;
    d_minimum = safe_d_multiplier * fzero(radiusMax_nn, 0);
    
    %Could additionally be defined as a maximum distance as with the upper
    %limit of the d-coefficient
   
    radiusMin_n = subs(radiusFunction - rMin, theta, theta_f/2);
    radiusMin_nn = @(d_max_in) double(subs(radiusMin_n, d, d_max_in));
    
    %Needs to be <1
    d_maximum = fzero(radiusMin_nn, [d_minimum, 1]);

    global theta_0 d_solution
    %solution = fmincon(@Tds_min, [(theta_0 + theta_f) * 0.5, d_solution], [], [], [], [], [], [], []);
%     solution = fmincon(@transferTimeOptimization, d_solution, [], [], [], [], [], [], []);
% 
%     value = transferTimeOptimization(solution);
%     
%     maximumTimeError = 1;
%     
%     global tf

%     if value > maximumTimeError
%         N = N + 1;
%         TOF_estimation = (1+N)*pi*sqrt((a_initial+a_final)^3/(8*mju));
%         tf = TOF_estimation;
%         tof_current = TOF_estimation;
%         
%         updateParameters();
%     end
end

function out_o = Tds_min(inputVector)
    %global thetaDotSquareFunction
    global timeFunction tof_current;
    syms theta d;

    theta_zero = inputVector(1);
    d_zero = inputVector(2);

    %out_o = abs(double(subs(subs(thetaDotSquareFunction, d, d_zero), theta, theta_zero)));
    out_o = abs(double(subs(subs(timeFunction, d, d_zero), theta, theta_zero)) - tof_current);

    fprintf("The minimum Error found: %e\n", out_o);
   
end
% 
% function minError_o = Te_min(inputVector);
%     global timeFunction
% 
%     syms theta d;
% 
% 
% end
% 
