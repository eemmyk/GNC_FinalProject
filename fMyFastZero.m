function [b,fval] = fMyFastZero(FunFcn, x, y, tols, paramVector, tof_current, theta_super)
                                               % Customised parameters are given to the function
tol = tols(1);
tolf = tols(2);


a = x(1);
b = x(2);

fa = y(1);%FunFcn(a, paramVector, tof_current, theta_super);
fb = y(2);%FunFcn(b, paramVector, tof_current, theta_super);

if ( fa == 0 )
    b = a;
    fval = fa;
    return
elseif ( fb == 0)
    % b = b;
    fval = fb;
    return
end

fc = fb;

% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end

    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end
        if p > 0
            q = -q;
        else
            p = -p;
        end
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
        else
            d = m;  e = m;
        end
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler
        b = b + d;
    elseif b > c
        b = b - toler;
    else
        b = b + toler;
    end

    fb = FunFcn(b, paramVector, tof_current, theta_super);

    %Additional TolFun exit check
    tolerf = 2.0*tolf*max(abs(b),1.0);
    if abs(fb) < tolerf
        break
    end
    
end % Main loop

fval = fb; % b is the best value
