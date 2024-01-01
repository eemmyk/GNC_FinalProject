function [deltaV_o] = optimalDVSolver(tof_in)
    global tof_current  theta_f currentTime plotDV_3D theta_vec;
    global d_solution d_minimum d_maximum initial_DeltaV;
    global timeFunction timeFunction_nn;

    global plotAccuracy theta1 theta2 theta_0;
    global thrustFunction thetaDotFunction;
    global radiusFunction thetaDotSquareFunction;
    
    tof_current = tof_in;
    
%     global N
%     N = 0;

    updateParameters(0)
%     safeTransferAngle = pi;

%     if theta_f < safeTransferAngle
%         deltaV_o = 1e24; %A big number
%         return;
%     end


%     timeFunction
%     d_minimum
%     d_maximum
%     theta1
%     theta2
%     theta_f
%     tof_current

    %Check if the tof is solvable within the limits of the d-coefficient
    t_error_min = trapz(theta_vec, timeFunction_nn(d_maximum, theta_vec)) - tof_current;
    t_error_max = trapz(theta_vec, timeFunction_nn(d_minimum, theta_vec)) - tof_current;
    
    if (t_error_min < 0) ~= (t_error_max < 0)
        opt = optimset('TolFun', 1e3);
        d_solution = fzero(@transferTimeOptimization, [d_minimum, d_maximum], opt);

        syms d theta;
        jerkFunction_nn = @(d_coeff, angle) double(subs(subs(thrustFunction/thetaDotFunction, d, d_coeff), theta, angle));

        deltaV_o = trapz(theta_vec, abs(jerkFunction_nn(d_solution, theta_vec)));
    else
        deltaV_o = 1e24; %A big number
        fprintf("Unfeasible Time Of Flight Detected\n")
    end

    global deltaResult 

    if deltaV_o < deltaResult
        %Get globals from updatedParameters
        global theta2 theta2_dot gamma2 r2 P2 theta1 theta1_dot gamma1 r1 P1
        global nu1_i nu2_i r1_i r2_i;
        
        %Save best results as globals
        global theta2_opt r2_opt tof_optimal theta2_dot_opt gamma2_opt P2_opt;
        global theta1_opt r1_opt theta1_dot_opt gamma1_opt P1_opt;
        global nu1_i_opt nu2_i_opt r1_i_opt r2_i_opt theta_f_opt d_opt;
        global thrustFunction_opt thetaDotFunction_opt thetaDotSquareFunction_opt;
        global radiusFunction_opt timeFunction_opt;

        deltaResult = deltaV_o;
        theta2_opt = theta2;
        theta1_opt = theta1;
        theta_f_opt = theta_f;
        theta2_dot_opt = theta2_dot;
        theta1_dot_opt = theta1_dot;
        gamma2_opt = gamma2;
        gamma1_opt = gamma1;
        d_opt = d_solution;
        tof_optimal = tof_in;
        r2_opt = r2;
        r1_opt = r1;
        P2_opt = P2;
        P1_opt = P1;
        nu1_i_opt = nu1_i;
        nu2_i_opt = nu2_i;
        r1_i_opt = r1_i;
        r2_i_opt = r2_i;

        thrustFunction_opt =  thrustFunction;
        thetaDotFunction_opt = thetaDotFunction; 
        thetaDotSquareFunction_opt = thetaDotSquareFunction; 
        timeFunction_opt = timeFunction;
        radiusFunction_opt = radiusFunction;


        fprintf("Optimal dV: %.0f m/s for transfer time: %.0f s\n", deltaV_o, tof_in)
    end

%     radiusFunction_n = subs(radiusFunction, d, d_solution);
%     radiusFunction_nn = @(angle) double(subs(radiusFunction_n, theta, angle));
%     theta_plot_vec = linspace(theta_0, theta_f, plotAccuracy);
% 
%     x = cos(theta_plot_vec+theta1) .* radiusFunction_nn(theta_plot_vec);
%     y = sin(theta_plot_vec+theta1) .* radiusFunction_nn(theta_plot_vec);
%    
%     plot(x, y, "Color", [0.2 0.7 0.2]);
   
    R_Multiplier = (3*(deltaV_o/initial_DeltaV)^2 - 2*(deltaV_o/initial_DeltaV)^3);
    G_Multiplier = 1-(3*((deltaV_o-initial_DeltaV)/initial_DeltaV)^2 - 2*((deltaV_o-initial_DeltaV)/initial_DeltaV)^3);


    if plotDV_3D == 1
        if deltaV_o > 2*initial_DeltaV
            color = [1 0 0];
        elseif deltaV_o > initial_DeltaV
            color = [1 G_Multiplier 0];
        elseif deltaV_o < initial_DeltaV
            color = [R_Multiplier 1 0];
        end
    %     
    % 
    %     G_Multiplier = 0.5 * (1 + cos(goodnessRatio * pi));
    %     R_Multiplier = 0.5 * (1 - cos(goodnessRatio * pi));
          
    %     color = [R_Multiplier G_Multiplier 0];

       
        plot3(currentTime, tof_current, deltaV_o,'-o','Color','b','MarkerSize',10,'MarkerFaceColor', color)
    else
        plot(tof_in, deltaV_o,'or', 'MarkerSize',2,'MarkerFaceColor','r')
    end
    pause(0.001);
end

