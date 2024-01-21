function [radiusVec] = fRadiusFunction(d, theta, paramVector)

    theta_f = paramVector.theta_f;
    r2 = paramVector.r2;
    
    a = paramVector.a;
    b = paramVector.b;
    c = paramVector.c;
                    
    efg_Mat_2 = [1/r2 - (a + b*theta_f + c*theta_f^2 + d*theta_f^3);
                -tan(paramVector.gamma2)/r2 - (b + 2*c*theta_f + 3*d*theta_f^2); 
                paramVector.mju/(r2^4*paramVector.theta2_dot^2) - (1/r2 + 2*c + 6*d*theta_f)];
    
    efg = 1/(2*theta_f^6) * paramVector.efg_Mat_1 * efg_Mat_2;
    
    e = efg(1);
    f = efg(2);
    g = efg(3);
    
    radiusVec = 1 ./ (a + b.*theta(1,:) + c.*theta(2,:) + d.*theta(3,:) + e.*theta(4,:) + f.*theta(5,:) + g.*theta(6,:));
end

