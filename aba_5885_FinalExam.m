function  aba_5885_FinalExam
    
    global m1 m2 L1 L2 K g delta_R delta_A t_0 t_f delta_t;
    constants = csvread('SU24FinalExamData.csv');
    m1 = constants(1);
    m2 = constants(2);
    L1 = constants(3);
    L2 = constants(4);
    K = constants(5);
    g = constants(6);
    delta_R = constants(7);
    delta_A = constants(8);
    t_0 = constants(9);
    t_f = constants(10);
    delta_t = constants(11);

    hFig = figure('Position', [100, 100, 800, 600], 'Name', 'Coupled Pendulum System');

    hAxes = axes('Parent', hFig, 'Units', 'normalized', 'Position', [0.1, 0.3, 0.6, 0.55]);
    title(hAxes, 'State Variables Responses Vs. Time');
    xlabel(hAxes, 'x axis label of Time (sec)');
    ylabel(hAxes, 'y axis label of State Variables Responses');
    xlim(hAxes, [-5 5]);
    ylim(hAxes, [-5 5]);
    grid(hAxes, 'on');

    hPanel = uipanel('Parent', hFig, 'Title', 'Initial Conditions', 'Position', [0.75, 0.4, 0.2, 0.5]);

    uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'theta1(0) in deg:', 'Position', [10, 270, 100, 20]);
    hTheta1 = uicontrol('Parent', hPanel, 'Style', 'edit', 'Position', [120, 270, 60, 20]);

    uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'theta1 dot (0) in deg/sec:', 'Position', [10, 220, 100, 20]);
    hTheta1Dot = uicontrol('Parent', hPanel, 'Style', 'edit', 'Position', [120, 220, 60, 20]);

    uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'theta2(0) in deg:', 'Position', [10, 170, 100, 20]);
    hTheta2 = uicontrol('Parent', hPanel, 'Style', 'edit', 'Position', [120, 170, 60, 20]);

    uicontrol('Parent', hPanel, 'Style', 'text', 'String', 'theta2 dot (0) in deg/sec:', 'Position', [10, 120, 100, 20]);
    hTheta2Dot = uicontrol('Parent', hPanel, 'Style', 'edit', 'Position', [120, 120, 60, 20]);

    % Plot and Reset buttons
    uicontrol('Parent', hFig, 'Style', 'pushbutton', 'String', 'Plot', 'Position', [1220, 200, 100, 40], ...
        'Callback', @(~,~) plotButtonCallback(hAxes, hTheta1, hTheta1Dot, hTheta2, hTheta2Dot));
    uicontrol('Parent', hFig, 'Style', 'pushbutton', 'String', 'Reset', 'Position', [1345, 200, 100, 40], ...
        'Callback', @(~,~) resetButtonCallback(hAxes));

    % Plot button function
    function plotButtonCallback(hAxes, hTheta1, hTheta1Dot, hTheta2, hTheta2Dot)

        theta1_0 = str2double(get(hTheta1, 'String'));
        theta1_dot_0 = str2double(get(hTheta1Dot, 'String'));
        theta2_0 = str2double(get(hTheta2, 'String'));
        theta2_dot_0 = str2double(get(hTheta2Dot, 'String'));

        % Angle limits
        if abs(theta1_0) > 11 || abs(theta2_0) > 11
            menu('Don''t forget the small angle approximation!', 'Okay');
            return;
        end

        % Deg to Radians
        theta1_0 = deg2rad(theta1_0);
        theta1_dot_0 = deg2rad(theta1_dot_0);
        theta2_0 = deg2rad(theta2_0);
        theta2_dot_0 = deg2rad(theta2_dot_0);

        % IC
        y0 = [theta1_0, theta1_dot_0, theta2_0, theta2_dot_0];

        options = odeset('RelTol', delta_R, 'AbsTol', delta_A, 'OutputFcn', @odeOutputFcn);

        tic;
        [t, y] = ode45(@pendulum_odes, t_0:delta_t:t_f, y0, options);
        elapsedTime = toc;

        plot(hAxes, t, y(:, 1), 'r-', t, y(:, 2), 'g--', t, y(:, 3), 'b:', t, y(:, 4), 'k-.');
        title(hAxes, 'State Variables Responses vs. Time');
        xlabel(hAxes, 'Time (sec)');
        ylabel(hAxes, 'State Variables Responses');
        legend(hAxes, '\theta_1', '\theta_1 dot', '\theta_2', '\theta_2 dot', 'Location', 'northeast');

        % Writing data to file
        dlmwrite('ode_solution.txt', [t, y], 'delimiter', ':');

        fid = fopen('ode_time.out', 'w');
        fprintf(fid, 'ODE45 elapsed time: %.4f seconds\n', elapsedTime);
        fclose(fid);
    end

    function resetButtonCallback(hAxes)
        cla(hAxes);
        xlim(hAxes, [-5 5]);
        ylim(hAxes, [-5 5]);
        title(hAxes, 'State Variables Responses vs. Time');
        xlabel(hAxes, 'Time (sec)');
        ylabel(hAxes, 'State Variables Responses');
        grid(hAxes, 'on');
    end


    function dydt = pendulum_odes(~, y)

        theta1 = y(1);
        theta1_dot = y(2);
        theta2 = y(3);
        theta2_dot = y(4);

        % Derivatives
        theta1_ddot = (theta1 * (m1 * (L1 * theta1_dot^2 - g) - K * L1) + K * L2 * theta2) / (m1 * L1);
        theta2_ddot = (theta2 * (m2 * (L2 * theta2_dot^2 - g) - K * L2) + K * L1 * theta1) / (m2 * L2);


        dydt = [theta1_dot; theta1_ddot; theta2_dot; theta2_ddot];
    end


    function status = odeOutputFcn(t, ~, flag)
        status = 0;
        if isempty(flag)
            if mod(t(end), delta_t) ~= 0
                status = 1;
            end
        end
    end
end