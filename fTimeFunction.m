function [timeStep] = fTimeFunction(d, theta, paramVector)

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

    timeStep = sqrt((r.^4./mju) .* (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4));
end