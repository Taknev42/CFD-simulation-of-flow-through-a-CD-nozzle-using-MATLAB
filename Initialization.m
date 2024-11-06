function [rho, T, A, v, P, throat_index] = Initialization(type, n, x)
    % Initialization of nozzle flow variables and area profile
    rho = zeros(1, n); % Density
    T = zeros(1, n);   % Temperature
    P = zeros(1, n);   % Pressure
    A = zeros(1, n);   % Area profile
    v = zeros(1, n);   % Velocity

    if type == "Subsonic_1"
        
        % Define the initial conditions along the nozzle
        for i = 1:n
            rho(i) = 1 - 0.023*x(i);
            T(i) = 1 - 0.009333*x(i);
            A(i) = 1 + 2.2 * (x(i) - 1.5)^2;  % Converging-diverging area profile
            v(i) = 0.59 / (rho(i) * A(i));      % Velocity
        end
        
        % Finding the index of the throat (Area is 1)
        throat_index = find(A == 1);
    elseif type == "Subsonic_2"
        
        
        % Define the initial conditions along the nozzle
        for i = 1:n
            % Linearly decreasing density and Temperature has been assumed
            % according to the initial conditions provided by John D. Anderson in
            % his book.
            rho(i) = 1 - 0.023*x(i);
            T(i) = 1 - 0.009333*x(i);
            if x(i) <= 1.5
                A(i) = 1 + 2.2 * (x(i) - 1.5)^2;  % Converging-diverging area profile
            elseif 1.5 <= x(i) <= 3
                A(i) = 1 + 0.2223*(x(i) - 1.5)^2;
            end
            v(i) = 0.59 / (rho(i) * A(i));      % Velocity
        
        end
        % plot(x, A)
        
        throat_index = find(A == 1);
    elseif type == "Supersonic"
        
        % Define the initial conditions along the nozzle. These conditions have
        % been picked from the ones defined by John D. Anderson in his book. 
        for i = 1:n
            if x(i) <= 0.5
                rho(i) = 1;       % Initial density before x = 0.5
                T(i) = 1;         % Initial temperature
            elseif (x(i) > 0.5 && x(i) <= 1.5)
                rho(i) = 1 - 0.366 * (x(i) - 0.5); % Linear decrease in density
                T(i) = 1 - 0.167 * (x(i) - 0.5);   % Linear decrease in temperature
            elseif (x(i) > 1.5 && x(i) <= 3)
                rho(i) = 0.634 - 0.3879 * (x(i) - 1.5); % Further decrease beyond x = 1.5
                T(i) = 0.833 - 0.3507 * (x(i) - 1.5);   % Further decrease in temperature
            end
        
            A(i) = 1 + 2.2 * (x(i) - 1.5).^2;  % Converging-diverging area profile
            v(i) = 0.59 / (rho(i) * A(i));      % Velocity
        end
        
        throat_index = find(A == 1);

    end
end
