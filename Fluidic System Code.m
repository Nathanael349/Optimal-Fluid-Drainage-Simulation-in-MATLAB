% Constants
g = 9.81;                               % Gravity (m/s^2)
rho = 1000;                             % Density of water (kg/m^3)
mu = 0.001;                             % Dynamic viscosity of water (Pa·s or kg/(m·s))
tank_height = 0.08;                     % Initial water height in the tank (m)
tank_area = 0.32 * 0.26;                % Cross-sectional area of tank (m^2)
tube_diameter = 0.00794;                % Tube diameter (m)
reservoir_diameter = 0.26;              % Reservoir diameter (m)
tube_radius = tube_diameter / 2;        % Tube radius (m)
tube_area = pi * tube_radius^2;         % Cross-sectional area of tube (m^2)
vertical_drop = 0.02;                   % Vertical drop for the tube exit point (m)
roughness = 0.0005;                     % 0.0025*10^-3         % Pipe roughness (m)


% Experimental tube lengths (in meters) and corresponding data in seconds
tube_lengths = [0.2, 0.3, 0.4, 0.6]; % Tube lengths (m)
experimental_times = [3*60+19, 3*60+34, 4*60+26, 4*60+48]; % Experimental drain times (s)

tube_test_quant = 15;

tube_lengths_test = zeros(1, tube_test_quant); % Preallocate array for efficiency
for i = 1:tube_test_quant
    tube_lengths_test(i) = 0.0 + (i-1)*0.1; % Start at 0.2, increment by 0.1 each step
end

num_lengths = tube_test_quant;

% Initialize arrays for results
drain_times = zeros(1, num_lengths);
distances = zeros(1, num_lengths);
angles = zeros(1, num_lengths); % Array to store calculated angles

% Loop over each tube length to calculate drain time and distance
for i = 1:num_lengths
    
    L = tube_lengths_test(i);
    
    % Initial conditions
    h = tank_height;     % Water height in the tank (m)
    dt = 0.01;          % Smaller time step for accuracy (s)
    t = 0;               % Initial time (s)

    while (h>0)
        v_exit_converged = find_v_exit(g,h,rho,mu,tube_diameter,roughness,L,reservoir_diameter);
        dH = (tube_area/tank_area)*v_exit_converged*dt;
        h = h - dH;
        t = t + dt;
    end

    % Record drain time
    drain_times(i) = t;

    % Calculate the angle theta based on the vertical drop and length L
    theta = asin(vertical_drop / L);
    angles(i) = rad2deg(theta); % Store theta in degrees for reference

    % Calculate horizontal distance using projectile motion
    distances(i) = v_exit_converged * cos(theta) * (2 * v_exit_converged * sin(theta) / g);
end

function [v_exit] = find_v_exit(g, tank_height, rho, mu, tube_diameter, roughness, tube_length, reservoir_diameter)

    tolerance = 1e-3;
    max_iterations = 100; % Set a limit for iterations
    iteration_count = 0;

    Re = 0;
    
    % initial guess [Torricelli] (+0.04 to account for exit point located 4cm below the end line)
    v_exit_original = sqrt(2 * g * (tank_height + 0.04));
    
    % initialize new velocity variable so its accessible from outside while loop
    v_exit_new = 0;

    while true

        % Reynolds number calculation
        Re = (rho * v_exit_original * tube_diameter) / mu;

        if Re > 4000
            % Fully turbulent (Haaland equation)
            f = 1 / (-1.8 * log10(((roughness / tube_diameter) / 3.7)^1.11 + 6.9 / Re))^2;
        elseif Re > 2300
            % Laminar friction factor
            f_laminar = 64 / Re;
            % Turbulent friction factor using Haaland equation
            f_turbulent = 1 / (-1.8 * log10(((roughness / tube_diameter) / 3.7)^1.11 + 6.9 / Re))^2;
            % Weight for laminar contribution
            w = (4000 - Re) / (4000 - 2300);
            % Weighted friction factor
            f = w * f_laminar + (1 - w) * f_turbulent;
        elseif Re > 0
            % Fully laminar
            f = 64 / Re;
        else
             error('Reynolds number is zero or negative. Check input parameters.');
        end

        % Calculate K value for minor loss (sudden contraction)
        K = 0.42 * (1 - tube_diameter^2 / reservoir_diameter^2); %0.42
        
        % Bernouilli (+0.04 to account for exit point located 4cm below the end line)
        % Incorporate head loss directly into Bernouilli
        v_exit_new = sqrt((2*g*(tank_height+0.04))/(1+f*tube_length/tube_diameter+K));

        % Check for convergence
        if abs(v_exit_new - v_exit_original) < tolerance
            break;
        end

        % Update for next iteration
        if iteration_count > max_iterations
            error('Velocity Convergence not achieved within maximum iterations');
        end

        v_exit_original = v_exit_new;
        iteration_count = iteration_count + 1;


    end
    
    % function output
    v_exit = v_exit_new;
     
end


% Display results in minutes:seconds format
disp('Tube Length (m)   Angle (degrees)   Drain Time (min:sec)   Horizontal Distance (m)')
for i = 1:num_lengths
    % Convert drain time to minutes and seconds
    minutes = floor(drain_times(i) / 60);
    seconds = round(mod(drain_times(i), 60));
    fprintf('%-15.2f %-18.2f %d:%02d %25.2f\n', tube_lengths_test(i), angles(i), minutes, seconds, distances(i));
end


% Indices for comparison
comparison_indices = [3, 4, 5, 7]; % Indices in drain_times that correspond to experimental_times

% Extract the relevant drain_times for comparison
model_times = drain_times(comparison_indices) / 60; % Convert to minutes

% Convert experimental_times to minutes
experimental_times_minutes = experimental_times / 60;

% Compute percentage errors
percentage_errors = abs(model_times - experimental_times_minutes) ./ experimental_times_minutes * 100;

% Find the maximum percentage error
max_percentage_error = max(percentage_errors);

% Compute the average percentage error
average_percentage_error = mean(percentage_errors);

fprintf('The maximum percentage error is %.2f%%.\n', max_percentage_error);
fprintf('The average percentage error is %.2f%%.\n', average_percentage_error);


% Plot results and compare with experimental data
figure; % Create a new figure
plot(tube_lengths_test, drain_times / 60, 'o-b', 'LineWidth', 2); % Convert time to minutes
hold on;
plot(tube_lengths, experimental_times / 60, 'x-r', 'LineWidth', 2); % Experimental data in minutes

% Add a vertical line at x = 0.5 (optimal value)
xline(0.5, '--k', 'LineWidth', 1.5); % Dashed black line with a width of 1.5

% Add labels, legend, and title
xlabel('Tube Length (m)');
ylabel('Drain Time (min)');
legend('Model Drain Time', 'Experimental Drain Time', 'Optimal');
title('Drain Time vs Tube Length');
hold off;



