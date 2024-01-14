function [timeStep] = fTimeFunction(d, theta, paramVector)

    theta1 = theta(1,:);
    theta2 = theta(2,:);
    theta3 = theta(3,:);
    theta4 = theta(4,:);
    theta5 = theta(5,:);
    theta6 = theta(6,:);

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
    
    r = 1 ./ (a + b.*theta1 + c.*theta2 + d.*theta3 + e.*theta4 + f.*theta5 + g.*theta6);

    timeStep = sqrt((r.^4./mju) .* (1./r + 2.*c + 6.*d.*theta1 + 12.*e.*theta2 + 20.*f.*theta3 + 30.*g.*theta4));
    
end