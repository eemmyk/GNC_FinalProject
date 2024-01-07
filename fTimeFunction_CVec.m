function [timeStep] = fTimeFunction_CVec(coeffVector, theta)
    global mju;
    a = coeffVector(1);
    b = coeffVector(2);
    c = coeffVector(3);
    d = coeffVector(4);
    e = coeffVector(5);
    f = coeffVector(6);
    g = coeffVector(7);
    
    r = 1 ./ (a + b.*theta + c.*theta.^2 + d.*theta.^3 + e.*theta.^4 + f.*theta.^5 + g.*theta.^6);
%     timeTroublePart = (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4);
% 
%     timeTroubleDiff = diff(timeTroublePart);
% 
%     timeTroubleZero = solve(timeTroubleDiff == 0, theta);
% 
%     timeTroubleRoot = vpa(timeTroubleZero)
% 

    timeStep = sqrt((r.^4./mju) .* (1./r + 2.*c + 6.*d.*theta + 12.*e.*theta.^2 + 20.*f.*theta.^3 + 30.*g.*theta.^4));

end
