function [time_t_f] = fTimeFunction(d, theta, dT, paramVector)

    theta_f = paramVector.theta_f;
    
    a = paramVector.a;
    b = paramVector.b;
    c = paramVector.c;

    efg_d_part = [-d*theta_f^3;
                  -3*d*theta_f^2;
                  -6*d*theta_f];

    efg = 1/(2*theta_f^6) * paramVector.efg_Mat_1 * (paramVector.efg_Mat_2_const + efg_d_part);
    
    e = efg(1);
    f = efg(2);
    g = efg(3);
    
    r = 1 ./ (a + b.*theta(1,:) + c.*theta(2,:) + d.*theta(3,:) + e.*theta(4,:) + f.*theta(5,:) + g.*theta(6,:));

    timeStep_Vec = (r.^2./sqrt(paramVector.mju)) .* sqrt((1./r + 2.*c + 6.*d.*theta(1,:) + 12.*e.*theta(2,:) + 20.*f.*theta(3,:) + 30.*g.*theta(4,:))); 

    time_t_f = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1));
end