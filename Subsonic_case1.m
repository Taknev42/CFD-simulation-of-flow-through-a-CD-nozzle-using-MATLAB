% Simulation of fluid flow through a converging-diverging nozzle using
% MacCormack method in conservation form. All equations are being solved in
% the non dimensionalised form.
clear all;
clc;
close all;

%% Numerical Solution
% Parameters
n = 200; % Number of grid points
x = linspace(0, 3, n); % Spatial domain (0 to 3)
dx = x(2) - x(1); % Grid spacing
g = 1.4; % Specific heat ratio (changed from gamma)
CFL = 0.5; % Courant number for stability
itr = 1000000; % Number of iterations (changed from n_iter)

% Initialization of primitive variables and Area profile
type = "Subsonic_1";
[rho, T, A, v, P, t_idx] = Initialization(type, n, x); % changed from throat_index
[Q1, Q2, Q3] = conserved(rho, A, v, T, g);
[F1, F2, F3] = flux(Q1, Q2, Q3, g);

%%
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
    
    % new Q obtained by averaging the derivatives
    Q1_n = Q1 + (0.5*dt)*(dQ1_dt_c + dQ1_dt_p); % changed from Q1_new
    Q2_n = Q2 + (0.5*dt)*(dQ2_dt_c + dQ2_dt_p); % changed from Q2_new
    Q3_n = Q3 + (0.5*dt)*(dQ3_dt_c + dQ3_dt_p); % changed from Q3_new
    
    % inlet boundary conditions
    Q1_n(1) = rho(1) * A(1);
    Q2_n(1) = 2*Q2_n(2) - Q2_n(3);
    v(1) = Q2_n(1)/Q1_n(1);
    Q3_n(1) = Q1_n(1) * (T(1) / (g - 1) + 0.5*g* v(1)^2);
    
    % outlet boundary conditions
    P(n) = 0.9999;
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
    
    % Extracting and Updating the primitive variables for the next iteration
    [rho, v, T, P, M] = update_primitives(Q1, Q2, Q3, A, g); % changed from Mach_no
    
    % Defining this condition so that we do not get more instabilities
    if M(t_idx) >= 0.8
        break
    end
end

% Plot results
figure;
plot(x, P, 'LineWidth', 1.5);
title('Pressure Distribution along the Nozzle');
xlabel('x/L');
ylabel('P/P_0');
grid on;

figure;
plot(x, T, 'LineWidth', 1.5);
title('Temperature Distribution along the Nozzle');
xlabel('x/L');
ylabel('T/T_0');
grid on;

figure;
plot(x, M, 'LineWidth', 1.5);
title('Mach Number Distribution along the Nozzle');
xlabel('x/L');
ylabel('M');
grid on;

figure;
plot(x, rho, 'LineWidth', 1.5);
title('Density Distribution along the Nozzle');
xlabel('x/L');
ylabel('\rho/\rho_0');
grid on;