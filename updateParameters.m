function [] = updateParameters()
    % Calculating the orbital parameters
    %initial orbit parameters
    global a_initial a_final currentTime mju N tof_current;
    global theta1 theta2 theta_f omega1 omega2 e1 e2 theta2_dot;
    global gamma1 r1 P1 gamma2 r2 P2 theta1_dot;
    global nu1_i nu2_i r1_i r2_i;



    a1 = a_initial;
    %calulating new point for object 1 as well
    %Time of last perigee pass for object 1
    Tp1 = 0;
    
    %Trying this one instead
    n1 = sqrt(mju/a1^3);
    P1 = 2*pi/n1;
    T1_i = mod(currentTime, P1) - Tp1;
    
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
    Tp2 = 0;
    
    %Trying this one instead
    n2 = sqrt(mju/a2^3);
    P2 = 2*pi/n2;
    T2_i = mod(currentTime, P2) - Tp2;
    
    syms nu_time
    
    E = 2*atan(tan((nu_time-pi)/2)/sqrt((1+e2)/(1-e2)));
    nuSolver = n2*T2_i==pi+E-e2*sin(E);
    nuSolutions_i = vpasolve(nuSolver, nu_time);
    nuSolutions_i = double(nuSolutions_i);
    
    nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);


    %Now we can continue by calculating the true anomaly when the spacecraft
    %reaches orbit 2. At currentTime + tof_current
    
    finalTime = currentTime + tof_current;
    T = mod(finalTime, P2) - Tp2;
    
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
    
    global thrustFunction thetaDotFunction thetaDotSquareFunction;

    thetaDotFunction = sqrt((mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4));
    thrustFunction = -mju / (2 * r^3 * cos(gamma)) * (6*d + 24*e*theta + 60*f*theta^2 + 120*g*theta^3 - tan(gamma)/r) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4)^2;
    
    thetaDotSquareFunction = (mju/r^4) / (1/r + 2*c + 6*d*theta + 12*e*theta^2 + 20*f*theta^3 + 30*g*theta^4);

end

