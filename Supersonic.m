% Simulation of fluid flow through a converging-diverging nozzle using MacCormack method
clear all;
clc;
close all;

%% Analytical Solution

% Parameters
g = 1.4;  % Specific heat ratio (changed from gamma)
x_a = linspace(0, 3, 1000);  % Spatial coordinates (changed from x_analytical)
A = 1 + 2.2*(x_a - 1.5).^2;  % Area distribution

% Find A* (throat area) which occurs at x = 1.5
A_t = 1;  % At x = 1.5, A = 1 (changed from A_star)

% Initialize arrays
M_a = zeros(size(x_a));  % Changed from M_analytical
P_r = zeros(size(x_a));  % Changed from P_P0
rho_r = zeros(size(x_a)); % Changed from rho_rho0
T_r = zeros(size(x_a));   % Changed from T_T0

% Function to solve for Mach number given area ratio (for supersonic flow)
func = @(M, A_ratio) A_ratio - (1./M).*(2/(g + 1)*(1 + (g-1)/2*M.^2)).^((g+1)/(2*(g-1)));

% Calculate Mach number for each point
for i = 1:length(x_a)
    A_ratio = A(i)/A_t;
    
    % Initial guess for Mach number
    if x_a(i) <= 1.5  % Subsonic region
        M_g = 0.5;    % Changed from M_guess
    else  % Supersonic region
        M_g = 2.0;
    end
    
    % Solve for Mach number
    opts = optimset('Display', 'off');  % Changed from options
    if x_a(i) <= 1.5  % Subsonic region
        M_a(i) = fsolve(@(M) func(M, A_ratio), M_g, opts);
    else  % Supersonic region
        M_a(i) = fsolve(@(M) func(M, A_ratio), M_g, opts);
    end
    
    % Calculating other properties
    P_r(i) = (1 + (g-1)/2*M_a(i)^2)^(-g/(g-1));
    rho_r(i) = (1 + (g-1)/2*M_a(i)^2)^(-1/(g-1));
    T_r(i) = (1 + (g-1)/2*M_a(i)^2)^(-1);
end


%% Numerical Solution

% Parameters
n = 300;              % Number of grid points
x = linspace(0, 3, n); % Spatial domain (0 to 3)
dx = x(2) - x(1);      % Grid spacing
g = 1.4;               % Specific heat ratio (gamma)
CFL = 0.4;             % Courant number for stability
itr = 10000;           % Number of time steps (changed from n_iter)

% Initialization of nozzle flow variables and area profile
type = "Supersonic";
[rho, T, A, v, P, t_idx] = Initialization(type, n, x);  % Changed throat_index to t_idx

% Initializing the conserved variables and Flux terms
[Q1, Q2, Q3] = conserved(rho, A, v, T, g);
[F1, F2, F3] = flux(Q1, Q2, Q3, g);

% Main time-stepping loop
for nt = 1:itr
    % Calculating time step using CFL condition
    dt_arr = zeros(1, n);
    for p = 1:n
        dt_arr(p) = (CFL*dx)/((T(p)^0.5)+v(p));
    end
    dt = min(dt_arr);
    check = isreal(dt);
    if (check ~= 1)
        break;
    end
    % predictor
    [Q1_p, Q2_p, Q3_p, dQ1_dt_p, dQ2_dt_p, dQ3_dt_p] = predictor_step(Q1, Q2, Q3, F1, F2, F3, dx, dt, A, rho, T, g);
    [F1_p, F2_p, F3_p] = flux(Q1_p, Q2_p, Q3_p, g);

    % corrector
    [Q1_c, Q2_c, Q3_c, dQ1_dt_c, dQ2_dt_c, dQ3_dt_c] = corrector_step(Q1_p, Q2_p, Q3_p, F1_p, F2_p, F3_p, dx, dt, A, rho, T, g);

    % new Q obtained by averaging the derivatives
    Q1_n = Q1 + (0.5*dt)*(dQ1_dt_c + dQ1_dt_p);  % Changed from Q1_new
    Q2_n = Q2 + (0.5*dt)*(dQ2_dt_c + dQ2_dt_p);  % Changed from Q2_new
    Q3_n = Q3 + (0.5*dt)*(dQ3_dt_c + dQ3_dt_p);  % Changed from Q3_new
    
    % Applying boundary conditions
    Q1_n(1) = rho(1) * A(1);
    Q2_n(1) = 2*Q2_n(2) - Q2_n(3);
    Q3_n(1) = Q1(1) * (T(1) / (g - 1) + 0.5*g* v(1)^2);
    
    % Extrapolating boundary conditions
    Q1_n(n) = 2*Q1_n(n-1) - Q1_n(n-2);
    Q2_n(n) = 2*Q2_n(n-1) - Q2_n(n-2);
    Q3_n(n) = 2*Q3_n(n-1) - Q3_n(n-1);
    
    % Updating Q values
    Q1 = Q1_n;
    Q2 = Q2_n;
    Q3 = Q3_n;

    [F1, F2, F3] = flux(Q1, Q2, Q3, g);
    
    % Extracting and Updating primitive variables
    [rho, v, T, P, M] = update_primitives(Q1, Q2, Q3, A, g);  % Changed Mach_no to M

end

% Plotting
figure;
plot(x_a, P_r, 'b-', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0 1 0.7]);
hold on;
plot(x, P, 'g--', 'LineWidth', 3, 'Color', [1 0 0 0.7]);
title('Pressure Distribution along the Nozzle');
xlabel('x/L');
ylabel('P/P_0');
legend('Analytical', 'Numerical', 'Location', 'best');
grid on;

figure;
plot(x_a, T_r, 'b-', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0 1 0.7]);
hold on;
plot(x, T, 'g--', 'LineWidth', 3, 'Color', [1 0 0 0.7]);
title('Temperature Distribution along the Nozzle');
xlabel('x/L');
ylabel('T/T_0');
legend('Analytical', 'Numerical', 'Location', 'best');
grid on;

figure;
plot(x_a, rho_r, 'b-', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0 1 0.7]);
hold on;
plot(x, rho, 'g--', 'LineWidth', 3, 'Color', [1 0 0 0.7]);
title('Density Distribution along the Nozzle');
xlabel('x/L');
ylabel('\rho/\rho_0');
legend('Analytical', 'Numerical', 'Location', 'best');
grid on;

figure;
plot(x_a, M_a, 'b-', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0 0 1 0.7]);
hold on;
plot(x, M, 'g--', 'LineWidth', 3, 'Color', [1 0 0 0.7]);
xlabel('x/L');
ylabel('M');
title('Mach Number Distribution');
legend('Analytical', 'Numerical', 'Location', 'best');
grid on;