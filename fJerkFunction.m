function [dV_o] = fJerkFunction(d, theta, dT, paramVector)

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
    gamma = atan(-r .* (b + 2.*c.*theta(1,:) + 3.*d.*theta(2,:) + 4.*e.*theta(3,:) + 5.*f.*theta(4,:) + 6.*g.*theta(5,:)));
    
    rsq = r.^2;

    thetaDotFunction = (sqrt(paramVector.mju)./rsq) ./ sqrt((1./r + 2.*c + 6.*d.*theta(1,:) + 12.*e.*theta(2,:) + 20.*f.*theta(3,:) + 30.*g.*theta(4,:)));
    thrustFunction = -paramVector.mju ./ (2.*rsq.*r.*cos(gamma)) .* (6.*d + 24.*e.*theta(1,:) + 60.*f.*theta(2,:) + 120.*g.*theta(3,:) - tan(gamma)./r) ./ (1./r + 2.*c + 6.*d.*theta(1,:) + 12.*e.*theta(2,:) + 20.*f.*theta(3,:) + 30.*g.*theta(4,:)).^2;
    
    acc_vec = abs(thrustFunction./thetaDotFunction);

    dV_o = dT * (acc_vec(1) + acc_vec(end)) / 2 + dT * sum(acc_vec(2:end-1));

end

