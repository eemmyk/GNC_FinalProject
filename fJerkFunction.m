function [jerk_vec] = fJerkFunction(d, theta)
    global mju r1 r2 gamma1 gamma2 theta_f theta1_dot theta2_dot
    
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
    gamma = atan(-r .* (b + 2.*c.*theta + 3.*d.*theta.^2 + 4.*e.*theta.^3 + 5.*f.*theta.^4 + 6.*g.*theta.^5));

    thetaDotFunction = sqrt((mju./r.^4) ./ (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4));
    thrustFunction = -mju ./ (2 .* r.^3 .* cos(gamma)) .* (6.*d + 24.*e.*theta + 60.*f.*theta.^2 + 120.*g.*theta.^3 - tan(gamma)./r) ./ (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4).^2;
    
    jerk_vec = thrustFunction./thetaDotFunction;

end

