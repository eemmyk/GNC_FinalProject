function [radiusVec] = fRadiusFunction(d, theta, useOptimal)
    global mju

    if useOptimal == 1
        global gamma1_opt gamma2_opt theta_f_opt 
        global theta1_dot_opt theta2_dot_opt r1_opt r2_opt;
        r1_l = r1_opt;
        r2_l = r2_opt;
        gamma1_l = gamma1_opt;
        gamma2_l = gamma2_opt;
        theta_f_l = theta_f_opt;
        theta1_dot_l = theta1_dot_opt;
        theta2_dot_l = theta2_dot_opt;
    else
        global gamma1 gamma2 theta_f theta1_dot theta2_dot r1 r2;
        r1_l = r1;
        r2_l = r2;
        gamma1_l = gamma1;
        gamma2_l = gamma2;
        theta_f_l = theta_f;
        theta1_dot_l = theta1_dot;
        theta2_dot_l = theta2_dot;
    end
    
    a = 1/r1_l;
    b = -tan(gamma1_l) / r1_l;
    c = 1/(2*r1_l) * (mju / (r1_l^3 * theta1_dot_l^2) - 1);
    
    efg_Mat_1 = [30*theta_f_l^2  -10*theta_f_l^3  theta_f_l^4;
                -48*theta_f_l     18*theta_f_l^2 -2*theta_f_l^3; 
                 20            -8*theta_f_l     theta_f_l^2];
    
    efg_Mat_2 = [1/r2_l - (a + b*theta_f_l + c*theta_f_l^2 + d*theta_f_l^3);
                -tan(gamma2_l)/r2_l - (b + 2*c*theta_f_l + 3*d*theta_f_l^2); 
                mju/(r2_l^4*theta2_dot_l^2) - (1/r2_l + 2*c + 6*d*theta_f_l)];
    
    efg = 1/(2*theta_f_l^6) * efg_Mat_1 * efg_Mat_2;
    
    e = efg(1);
    f = efg(2);
    g = efg(3);
    
    radiusVec = 1 ./ (a + b.*theta + c.*theta.^2 + d.*theta.^3 + e.*theta.^4 + f.*theta.^5 + g.*theta.^6);
end
