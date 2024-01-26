%Define starting true anomaly
%nu2_i = 0;

%Time of last perigee pass for object 2
Tp2 = 0;
currentTime = 90000;

%This can be used to calculate the current true anomaly with:
% M = nT = E - e*sin(E)
% sin(E) = sin(nu2_tf)*sqrt(1-e2^2) / (1+e2*cos(nu2_tf))
% n2 = sqrt(mju/a2^3);
% T = currentTime - Tp2;
% 
% syms nu2_i
% 
% nuSolver_i = n2*T==asin(sin(nu2_i)*sqrt(1-e2^2)/(1+e2*cos(nu2_i)))-e2*sin(nu2_i)*sqrt(1-e2^2) / (1+e2*cos(nu2_i));
% nuSolutions_i = vpasolve(nuSolver_i, nu2_i);
% nuSolutions_i = double(nuSolutions_i);



%Trying this one instead
n2 = sqrt(mju/a2^3);
P2 = 2*pi/n2;
T = mod(currentTime, P2) - Tp2;

% M = nT = E - e*sin(E)
% E = 2 * atan(tan((nu-pi)/2)/sqrt((1+e2)/(1-e2)));

% n2 = sqrt(mju/a2^3);
% T = currentTime - Tp2;


syms nu2_i

E = 2*atan(tan((nu2_i-pi)/2)/sqrt((1+e2)/(1-e2)));
nuSolver_i = n2*T==pi+E-e2*sin(E);
nuSolutions_i = vpasolve(nuSolver_i, nu2_i);
nuSolutions_i = double(nuSolutions_i);

nu2_i = mod(nuSolutions_i + 2*pi, 2*pi);

%%
%Solve nu at tf based on kelper's equations:
%Calculating the true anomaly at t = tf
% (M-M0)/n = t-t0
% n = sqrt(mju / a^3)
% M = E-e*sin(E)
% E = acos((e+cos(nu)) / (1+e*cos(nu)))

% sinh(H) = sin(nu) * sqrt(e2^2-1) / (1 + e*cos(nu))

% Xi = e2*sin(nu2_i)*sqrt(e2^2-1)/(1+e2*cos(nu2_i))-asinh(sin(nu2_i)*sqrt(e2^2-1)/(1 + e2*cos(nu2_i)));
% 
% syms nu2_tf
% 
% nuSolver = sqrt(-a2^3 /mju)*(e2*sin(nu2_tf)*sqrt(e2^2-1)/(1 + e2*cos(nu2_tf))-asinh(sin(nu2_tf)*sqrt(e2^2-1)/(1 + e2*cos(nu2_tf)))-Xi) == tf;
% nuSolutions = solve(nuSolver, nu2_tf);
% nuSolutions_n = double(nuSolutions);
% 
% n2 = sqrt(mju/a2^3);
% M0 = acos((e2+cos(nu2_i))/(1+e2*cos(nu2_i)))-e2*(sin(nu2_i)*sqrt(1-e2^2)/(1+e2*cos(nu2_i)));
% 
% syms nu2_tf
% nuSolver = ((acos((e2+cos(nu2_tf))/(1+e2*cos(nu2_tf)))-e2*(sin(nu2_tf)*sqrt(1-e2^2)/(1+e2*cos(nu2_tf))))-M0)/n2 == tf;
% nuSolutions = vpasolve(nuSolver, nu2_tf);
% nuSolutions_n = double(nuSolutions);
% %Since no retrograde orbits are considered, select the positive nu
% nu2 = nuSolutions_n > 0
% 
% % M = 2*pi*(t)

% The radius of the orbit can be represented as:
% r(theta) 1/(a + b*theta + c*theta.^2 + d*theta.^3 + e*theta.^4 + f*theta.^5)

%theta = theta_0 = 0 yields:
%r1 = 1/a

%theta = theta_f yields:
%r2 = 1/(a + b*theta_f + c*theta_f^2 + d*theta_f^3 + e*theta_f^4 + f*theta_f^5);

% The time derivative of the radius equation yields
% r_dot = -r.^2 * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)*theta_dot
% From this we can get the flight-path angle gamma may be found
% tan(gamma) = r_dot / (r * theta_dot) = -r * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)


%The first and final flight path angles gamma1 and gamma2 can be solved from this equation:
%tan(gamma1) = -r1 * b
%tan(gamma2) = -r2 * (b + 2*c*theta + 3*d*theta.^2 + 4*e*theta.^3 + 5*f*theta.^4)

%Additionally the polynomial must satisfy the equations of motion
%r_dot_dot = r*theta_dot^2 + mju/r^2 = Ta * sin(alpha)
%1 / r d/dt (r^2*theta_dot) = Ta * cos(alpha)

%Ta is the thrust acceleration and alpha is the thrust angle

%r*theta_dot_dot + 2*r*tan(gamma) = Ta * cos(alpha)

%A LOT OF STUFF HERE:


%%%%%%%%%%%%%%%%%%%%



function [deltaV_o] = optimalDVSolver(inputVec, pSettings)
    global d_solution deltaResult theta_super pState; 

    if pSettings.solveDate == 1
        pState.currentTime = inputVec(1);
        pState.tof_current = inputVec(2);
    else
        pState.tof_current = inputVec;
    end

    trueSolution = 0;
    N_start = pState.N;
    localBestDV = Inf;

    %Set limits for multiorbit search
    if pSettings.useMultiorbitFilling == 1
        N_end = pSettings.maxDepthN;
    else
        N_end = N_start;
    end

    pState.testedOrbits = pState.testedOrbits + 1;
    
    d_solution = 0;
    deltaV_o = 1e24; %A big number

    for N_current = N_start:N_end
        
        pState.N = N_current;
        [resultVector, paramVector] = updateParameters(0, pSettings);
        dT = theta_super(1,2) - theta_super(1,1);
        tof_current = pState.tof_current;

        d_minimum = resultVector(1);
        d_maximum = resultVector(2);
        realOrbit = resultVector(4);

%         figure(7);
%         hold on;
%         

%         for d_i = linspace(d_minimum, d_maximum, 1000)
%             timeStep_Vec = fTimeFunction(d_i, theta_super, paramVector);
%             f_i = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
%             plot(d_i, f_i, 'o')
%         end


     if ~realOrbit
        continue
     end

            x_min = d_minimum;
            timeStep_Vec = fTimeFunction(x_min, theta_super, paramVector);
            f_min = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
            x_max = d_maximum;
            timeStep_Vec = fTimeFunction(x_max, theta_super, paramVector);
            f_max = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
            
            crossing = (f_min < 0) ~= (f_max < 0);
            if ~crossing % || imaginary
                continue
            end

         try
            tfTimeHandle = @(d_in) transferTimeSolution(d_in, paramVector, pState.tof_current, theta_super);
            d_solution = fMyFastZero(tfTimeHandle, [d_minimum, d_maximum], pSettings.opt_tof_fzero);

%             timeStep_Vec = fTimeFunction(d_minimum, theta_super, paramVector);
%             f_min = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% 
%             timeStep_Vec = fTimeFunction(d_maximum, theta_super, paramVector);
%             f_max = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current; 


%             tfTimeRootHandle = @(d_in) fTimeFunction(d_in, theta_super, paramVector);
%             d_solution = myFzero(tfTimeRootHandle, d_minimum, d_maximum);
% 


            dT = theta_super(1,2) - theta_super(1,1);
            deltaV_o_Vec = abs(fJerkFunction(d_solution, theta_super, paramVector));
            deltaV_o = dT * (deltaV_o_Vec(1) + deltaV_o_Vec(end)) / 2 + dT * sum(deltaV_o_Vec(2:end-1));
            trueSolution = 1;
%         end
%             if realOrbit
%     
%                 x_min = d_minimum;
%                 timeStep_Vec = fTimeFunction(x_min, theta_super, paramVector);
%                 f_min = real(dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1))) - tof_current;
%                 x_max = d_maximum;
%                 timeStep_Vec = fTimeFunction(x_max, theta_super, paramVector);
%                 f_max = real(dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1))) - tof_current;
%                 
%                 crossing = (f_min < 0) ~= (f_max < 0);
%                 %imaginary = ~isreal(f_min) || ~isreal(f_max);
% 
%                 if ~crossing % || imaginary
%                     continue
%                 end
% 
% % 
% %                 if (d_minimum < 0) ~= (d_maximum < 0)
% %                     timeStep_Vec = fTimeFunction(0, theta_super, paramVector);
% %                     f_zero = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% %                     if f_zero > 0
% %                         f_min = f_zero;
% %                         x_min = 0;
% %                     else
% %                         f_max = f_zero;
% %                         x_max = 0;
% %                     end
% %                 end
% % 
% %                 %predict root
% %                 x1 = 0.99 * x_min + 0.01 * x_max;
% %                 x2 = 0.01 * x_min + 0.99 * x_max;
% % 
% %                 timeStep_Vec = fTimeFunction(x1, theta_super, paramVector);
% %                 y1 = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% % 
% %                 timeStep_Vec = fTimeFunction(x2, theta_super, paramVector);
% %                 y2 = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% %                 
% %                 k = (x1-d_minimum) / (x2-d_minimum);
% % 
% %                 b = (y1 * k + y2) / (1 - k);
% %                 a = (y1 - b)*(x1-d_minimum);
% % 
% %                 x_guess = (a-b*d_minimum) / -b;
% % 
% %                 timeStep_Vec = fTimeFunction(x_guess, theta_super, paramVector);
% %                 y_guess = real(dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1))) - tof_current;
% % 
% %                 
% %                 if abs(f_min) < abs(f_max)
% %                     fnm1 = f_min;
% %                     xnm1 = x_min;
% %                 else
% %                     fnm1 = f_max;
% %                     xnm1 = x_max;
% %                 end
% % 
% %                 fn = y_guess;
% %                 xn = x_guess;
% 
% %                 iter = 0;
% % 
% %                 %Bisection
% %                 while iter < 6
% %                     iter = iter + 1;
% % 
% %                     x_new = (x_min + x_max)/2;
% %                     %x_new = (x_min * f_max - x_max * f_min) / (f_min - f_max);
% %                     timeStep_Vec = fTimeFunction(x_new, theta_super, paramVector);
% %                     f_new = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% %     
% %                     if f_new > 0
% % %                         xnm1 = x_min;
% % %                         xn = x_new;
% % %                         fnm1 = f_min;
% % %                         fn = f_new;
% % 
% %                         x_min = x_new;
% %                         f_min = f_new;
% %                     else
% % %                         xnm1 = x_max;
% % %                         xn = x_new;
% % %                         fnm1 = f_max;
% % %                         fn = f_new;
% %             
% %                         x_max = x_new;
% %                         f_max = f_new;
% %                     end
% %                 end
% 
% %                 xnm1 = x_min;
% %                 xn = x_max;
% %                 fnm1 = f_min;
% %                 fn = f_max;
%     
%                 iter = 0;
% 
%                 f_new = 2;
%     
%                 while abs(f_new) > 1
%                     x_new = x_max - f_max * (x_max - x_min) / (f_max - f_min);
% 
% %                         if xnp1 < d_minimum
% %                             xnp1 = d_minimum;
% %                         elseif xnp1 > d_maximum
% %                             xnp1 = d_maximum;
% %                         end
% 
%                     timeStep_Vec = fTimeFunction(x_new, theta_super, paramVector);
%                     f_new = dT * (timeStep_Vec(1) + timeStep_Vec(end)) / 2 + dT * sum(timeStep_Vec(2:end-1)) - tof_current;
% 
%                     if f_new > 0
%                         x_min = x_new;
%                         f_min = f_new;
%                     else
%                         x_max = x_new;
%                         f_max = f_new;
%                     end
%                         
%                     iter = iter + 1;
%     
%                     if iter > 20 %|| ~isreal(fn)
%                         break;
%                     end
%                 end
%     
%                 if ~isreal(f_new)
%                     continue
%                 end
%     
%                 d_solution = x_new;
%                 deltaV_o_Vec = abs(fJerkFunction(d_solution, theta_super, paramVector));
%                 deltaV_o = dT * (deltaV_o_Vec(1) + deltaV_o_Vec(end)) / 2 + dT * sum(deltaV_o_Vec(2:end-1));
%     
%                 trueSolution = 1;
% % 
% % %                 if iter < 50
% % %                     d_solution = x_new;
% % %     %             
% % %                     %dT = theta_super(1,2) - theta_super(1,1);
% % %                     deltaV_o_Vec = abs(fJerkFunction(d_solution, theta_super, paramVector));
% % %                     deltaV_o = dT * (deltaV_o_Vec(1) + deltaV_o_Vec(end)) / 2 + dT * sum(deltaV_o_Vec(2:end-1));
% % %         
% % %                     %deltaV_o = trapz(theta_super(1,:), abs(fJerkFunction(d_solution, theta_super, paramVector)));
% % %                     trueSolution = 1;
% % %                 else
% % %                     d_solution = 0;
% % %                     deltaV_o = 1e24; %A big number
% % %                 end
%             else
%                 d_solution = 0;
%                 deltaV_o = 1e24; %A big number
%             end      
        %else
        catch
                d_solution = 0;
                deltaV_o = 1e24; %A big number
        end

        if deltaV_o <= localBestDV
            localBestDV = deltaV_o;
        %else
            %break
        end

        if deltaV_o < deltaResult
            %Save best results as globals
            global dateOptimal tof_optimal d_opt;
            global paramVector_opt
    
            deltaResult = deltaV_o;
            d_opt = d_solution;
    
            dateOptimal = pState.currentTime;
            tof_optimal = pState.tof_current;
    
            paramVector_opt = paramVector;
        end

    end

    if pSettings.plotTransferWindow == 1
        colorScaleDown = 4;

        R_Multiplier = (3*(localBestDV/pState.initial_DeltaV)^2 - 2*(localBestDV/pState.initial_DeltaV)^3);
        G_Multiplier = 1-(3*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^2 - 2*((localBestDV-pState.initial_DeltaV)/(pState.initial_DeltaV * colorScaleDown))^3);
    
        if localBestDV > (1+colorScaleDown)*pState.initial_DeltaV
            color = [1 0 0];
        elseif localBestDV > pState.initial_DeltaV
            color = [1 G_Multiplier^4 0];
        elseif localBestDV <= pState.initial_DeltaV
            color = [R_Multiplier^4 1 0];
        end

        color = color.*(0.75 + 0.25 * trueSolution);

        rectangle('Position',[pState.currentTime-0.5*pSettings.tfWindowPixelsX, pState.tof_current-0.5*pSettings.tfWindowPixelsY, ...
                              pSettings.tfWindowPixelsX, pSettings.tfWindowPixelsY], 'FaceColor', color, 'EdgeColor', color);
                           
    end   

    %Reset value
    pState.N = N_start;
    
    if trueSolution == 0
        pState.failedOrbits = pState.failedOrbits + 1;
    end
end

% function root = myFzero(func, x0, x1, tol, max_iter)
%     global theta_super pState
% 
%     dT = theta_super(1,2) - theta_super(1,1);
% 
% 
%     if nargin < 4
%         tol = 1e-6;
%     end
%     
%     if nargin < 5
%         max_iter = 100;
%     end
%     
%     for iteration = 1:max_iter
%         f0_vec = func(x0);
%         f0 = dT * (f0_vec(1) + f0_vec(end)) / 2 + dT * sum(f0_vec(2:end-1)) - pState.tof_current;
% 
%         f1_vec = func(x1);
%         f1 = dT * (f1_vec(1) + f1_vec(end)) / 2 + dT * sum(f1_vec(2:end-1)) - pState.tof_current;
%         
%         if abs(f1) < tol
%             root = x1;
%             return;
%         end
%         
%         % Check for sign change and apply bisection
%         if (f0 < 0) ~= (f1 < 0)
%             x2 = (x0 + x1) / 2;
%         else
%             % Apply secant method
%             x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
%             
%             % Check for convergence
% %             if abs(x2 - x1) < tol
% %                 root = x2;
% %                 return;
% %             end
%             
%             % Apply inverse quadratic interpolation
%             x3 = x0 - (f0 * (x1 - x0)^2) / (f1 - 2*f0 + f1);
% %             if abs(x3 - x2) < tol
% %                 root = x3;
% %                 return;
% %             end
% 
%             f3_vec = func(x3);
%             f3 = dT * (f3_vec(1) + f3_vec(end)) / 2 + dT * sum(f3_vec(2:end-1)) - pState.tof_current;
%             
%             % Choose the best approximation based on the bisection method
%             if (f1 < 0) ~= (f3 < 0)
%                 x2 = (x1 + x3) / 2;
%             end
%         end
%         
%         % Update iterates
%         x0 = x1;
%         x1 = x2;
%     end
%     
%     error('Root not found within the maximum number of iterations.');
% end