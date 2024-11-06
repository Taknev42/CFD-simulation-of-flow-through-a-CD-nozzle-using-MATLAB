% Simulation of fluid flow through a converging-diverging nozzle using MacCormack method
clear all;
clc;
close all;

%% Analytical Solution

g = 1.4;           % Specific heat ratio (changed from gamma)
P_o = 0.93;        % Outlet pressure ratio (changed from P_out)
nx = 100;          % Number of grid points for analytical solution (changed from nx_analytical)
x_a = linspace(0, 3, nx); % Spatial domain (changed from x_analytical)

% Area profile for analytical solution
A_a = zeros(1, nx);   % changed from A_analytical
for i = 1:nx
    if x_a(i) <= 1.5
        A_a(i) = 1 + 2.2 * (x_a(i) - 1.5)^2;
    else
        A_a(i) = 1 + 0.2223 * (x_a(i) - 1.5)^2;
    end
end

% Calculate outlet Mach number from pressure ratio
M_o = sqrt((2/(g-1))*((P_o)^(-(g-1)/g) - 1));  % changed from M_out
A_o = A_a(end);   % changed from A_out
A_t = A_o/((1/M_o) * ((2+(g-1)*M_o^2)/(g+1))^((g+1)/(2*(g-1)))); % changed from A_star

% Initialize arrays for analytical solution
M_a = zeros(1, nx);    % changed from M_analytical
P_a = zeros(1, nx);    % changed from P_analytical
T_a = zeros(1, nx);    % changed from T_analytical
rho_a = zeros(1, nx);  % changed from rho_analytical

% Calculate analytical solution
for i = 1:nx
    A_r = A_a(i)/A_t;    % changed from A_ratio
    M_g = 0.1;           % changed from M_guess
    
    % Newton-Raphson iteration
    M_c = M_g;           % changed from M_current
    err = 1;             % changed from error
    max_i = 1000;        % changed from max_iter
    nt = 0;
    
    while err > 1e-6 && nt < max_i
        fn = (1/M_c) * ((2+(g-1)*M_c^2)/(g+1))^((g+1)/(2*(g-1))) - A_r;     % changed from func
        dfn = -1/M_c^2 * ((2+(g-1)*M_c^2)/(g+1))^((g+1)/(2*(g-1))) + ...    % changed from dfunc
              (1/M_c) * ((2+(g-1)*M_c^2)/(g+1))^((g+1)/(2*(g-1))) * ...
              ((g+1)/(2*(g-1))) * (2*(g-1))/(2+(g-1)*M_c^2) * 2*M_c;
        M_n = M_c - fn/dfn;    % changed from M_new
        err = abs(M_n - M_c);
        M_c = M_n;
        nt = nt + 1;
    end
    
    M_a(i) = M_c;
end

% Calculate other analytical properties
P_a = (1 + (g-1)/2 * M_a.^2).^(-g/(g-1));
T_a = (1 + (g-1)/2 * M_a.^2).^(-1);
rho_a = (1 + (g-1)/2 * M_a.^2).^(-1/(g-1));

%% Numerical Solution

% Parameters
n = 200;              % Number of grid points
x = linspace(0, 3, n); % Spatial domain
dx = x(2) - x(1);      % Grid spacing
g = 1.4;               % Specific heat ratio
CFL = 0.5;             % Courant number for stability
itr = 100000;          % Number of iterations (changed from n_iter)

type = "Subsonic_2";
[rho, T, A, v, P, t_idx] = Initialization(type, n, x);  % changed from throat_index

[Q1, Q2, Q3] = conserved(rho, A, v, T, g);
[F1, F2, F3] = flux(Q1, Q2, Q3, g);

for nt = 1:itr
    % dt
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

    % new Q
    Q1_n = Q1 + (0.5*dt)*(dQ1_dt_c + dQ1_dt_p);   % changed from Q1_new
    Q2_n = Q2 + (0.5*dt)*(dQ2_dt_c + dQ2_dt_p);   % changed from Q2_new
    Q3_n = Q3 + (0.5*dt)*(dQ3_dt_c + dQ3_dt_p);   % changed from Q3_new

    % inlet boundary conditions
    Q1_n(1) = rho(1) * A(1);
    Q2_n(1) = 2*Q2_n(2) - Q2_n(3);
    v(1) = Q2_n(1)/Q1_n(1);
    Q3_n(1) = Q1_n(1) * (T(1) / (g - 1) + 0.5*g* v(1)^2);

    % outlet boundary conditions
    P(n) = 0.93;
    Q1_n(n) = 2*Q1_n(n-1) - Q1_n(n-2);
    Q2_n(n) = 2*Q2_n(n-1) - Q2_n(n-2);
    v(n) = Q2_n(n)/Q1_n(n);
    Q3_n(n) = P(n)*A(n)/(g - 1) + 0.5*g*Q2_n(n)*v(n);

    % Updating Q
    Q1 = Q1_n;
    Q2 = Q2_n;
    Q3 = Q3_n;

    % new F
    [F1, F2, F3] = flux(Q1, Q2, Q3, g);

    % Extracting new primitive variables for the next iteration
    [rho, v, T, P, M] = update_primitives(Q1, Q2, Q3, A, g);   % changed from Mach_no
end

figure;
plot(x_a, M_a, 'b-', 'LineWidth', 1.5);
hold on;
plot(x, M, 'r--', 'LineWidth', 1.5);
xlabel('x/L');
ylabel('Mach Number');
title('Mach Number Distribution');
legend('Analytical', 'Numerical', 'Location', 'best');
grid on;

% Area profile on second y-axis
yyaxis right;
plot(x_a, A_a, 'k:', 'LineWidth', 1);
ylabel('Area Ratio (A/A_{ref})');
legend('Analytical', 'Numerical', 'Area Profile', 'Location', 'best');

% Figure 2: Flow Properties Comparison
figure;
plot(x_a, P_a, 'b-', x_a, T_a, 'g-', x_a, rho_a, 'm-', 'LineWidth', 1.5);
hold on;
plot(x, P, 'b--', x, T, 'g--', x, rho, 'm--', 'LineWidth', 1.5);
xlabel('x/L');
ylabel('Property Ratio');
title('Flow Properties Distribution');
legend('P/P_0 (Analytical)', 'T/T_0 (Analytical)', '\rho/\rho_0 (Analytical)', ...
       'P/P_0 (Numerical)', 'T/T_0 (Numerical)', '\rho/\rho_0 (Numerical)', ...
       'Location', 'best');
grid on;