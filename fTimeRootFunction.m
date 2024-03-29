function [timeStep, timeDerivativeStep] = fTimeRootFunction(d, theta, paramVector)

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
    
    timeDerivativeStep = ((6.*theta1 + theta3 - (36.*theta2)./theta_f + (60.*theta3)./theta_f.^2 - (3.*theta4)./theta_f - (30.*theta4)./theta_f.^3 + (3.*theta5)./theta_f.^2 - theta6./theta_f.^3)./(mju.*(d.*theta3 + 1./r1 - theta6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) - theta4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) + theta5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) + (theta2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1) - (theta1.*tan(gamma1))./r1).^4) - (4.*(theta3 - (3.*theta4)./theta_f + (3.*theta5)./theta_f.^2 - theta6./theta_f.^3).*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta1 + d.*theta3 - theta4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) + (15.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (5.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) - theta6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) + (10.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (4.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + theta5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 + (24.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (9.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) - theta2.*((6.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^2 + (180.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^4 - (60.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^3) - theta4.*((15.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^4 + (300.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^6 - (120.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^5) + theta3.*((20.*((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2)))./theta_f.^3 + (480.*(d.*theta_f.^3 + 1./r1 - 1./r2 - (theta_f.*tan(gamma1))./r1 + (theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1)))./theta_f.^5 - (180.*(3.*d.*theta_f.^2 + ((mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - tan(gamma1)./r1 + tan(gamma2)./r2))./theta_f.^4) + 1./r1 + (theta2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1) - (theta1.*tan(gamma1))./r1))./(mju.*(d.*theta3 + 1./r1 - theta6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) - theta4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) + theta5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) + (theta2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1) - (theta1.*tan(gamma1))./r1).^5))./(2.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta1 + d.*theta3 - theta2.*(((6.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 36.*d.*theta_f + 6./r2 - (6.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^2 - (180.*d.*theta_f.^2 + (60.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (60.*tan(gamma1))./r1 + (60.*tan(gamma2))./r2)./theta_f.^3 + (180.*d.*theta_f.^3 + (90.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (180.*tan(gamma1).*theta_f)./r1 + 180./r1 - 180./r2)./theta_f.^4) - theta4.*(((15.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 90.*d.*theta_f + 15./r2 - (15.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^4 - (360.*d.*theta_f.^2 + (120.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (120.*tan(gamma1))./r1 + (120.*tan(gamma2))./r2)./theta_f.^5 + (300.*d.*theta_f.^3 + (150.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (300.*tan(gamma1).*theta_f)./r1 + 300./r1 - 300./r2)./theta_f.^6) + theta3.*(((20.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1 + 120.*d.*theta_f + 20./r2 - (20.*mju)./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (540.*d.*theta_f.^2 + (180.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (180.*tan(gamma1))./r1 + (180.*tan(gamma2))./r2)./theta_f.^4 + (480.*d.*theta_f.^3 + (240.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (480.*tan(gamma1).*theta_f)./r1 + 480./r1 - 480./r2)./theta_f.^5) + 1./r1 - theta6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - (12.*d.*theta_f.^2 + (4.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (4.*tan(gamma1))./r1 + (4.*tan(gamma2))./r2)./theta_f.^5 + (10.*d.*theta_f.^3 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (10.*tan(gamma1).*theta_f)./r1 + 10./r1 - 10./r2)./theta_f.^6) - theta4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - (15.*d.*theta_f.^2 + (5.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (5.*tan(gamma1))./r1 + (5.*tan(gamma2))./r2)./theta_f.^3 + (15.*d.*theta_f.^3 + (15.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./(2.*r1) - (15.*tan(gamma1).*theta_f)./r1 + 15./r1 - 15./r2)./theta_f.^4) + theta5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - (27.*d.*theta_f.^2 + (9.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f)./r1 - (9.*tan(gamma1))./r1 + (9.*tan(gamma2))./r2)./theta_f.^4 + (24.*d.*theta_f.^3 + (12.*(mju./(r1.^3.*theta1_dot.^2) - 1).*theta_f.^2)./r1 - (24.*tan(gamma1).*theta_f)./r1 + 24./r1 - 24./r2)./theta_f.^5) + (theta2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1) - (theta1.*tan(gamma1))./r1)./(mju.*(d.*theta3 + 1./r1 - theta6.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^4) - ((4.*tan(gamma2))./r2 - (4.*tan(gamma1))./r1 + 12.*d.*theta_f.^2 + (4.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^5 + (10.*d.*theta_f.^3 + 10./r1 - 10./r2 - (10.*theta_f.*tan(gamma1))./r1 + (5.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^6) - theta4.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./(2.*theta_f.^2) - ((5.*tan(gamma2))./r2 - (5.*tan(gamma1))./r1 + 15.*d.*theta_f.^2 + (5.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^3 + (15.*d.*theta_f.^3 + 15./r1 - 15./r2 - (15.*theta_f.*tan(gamma1))./r1 + (15.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1))./theta_f.^4) + theta5.*(((mju./(r1.^3.*theta1_dot.^2) - 1)./r1 + 6.*d.*theta_f + 1./r2 - mju./(r2.^4.*theta2_dot.^2))./theta_f.^3 - ((9.*tan(gamma2))./r2 - (9.*tan(gamma1))./r1 + 27.*d.*theta_f.^2 + (9.*theta_f.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^4 + (24.*d.*theta_f.^3 + 24./r1 - 24./r2 - (24.*theta_f.*tan(gamma1))./r1 + (12.*theta_f.^2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./r1)./theta_f.^5) + (theta2.*(mju./(r1.^3.*theta1_dot.^2) - 1))./(2.*r1) - (theta1.*tan(gamma1))./r1).^4)).^(1./2)); 

end