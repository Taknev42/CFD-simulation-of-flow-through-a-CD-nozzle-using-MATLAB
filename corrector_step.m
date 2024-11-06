function [Q1_c, Q2_c, Q3_c, dQ1_dt_c, dQ2_dt_c, dQ3_dt_c] = corrector_step(Q1_p, Q2_p, Q3_p, F1_p, F2_p, F3_p, dx, dt, A, rho, T, gamma)
    n = length(Q1_p);
    dQ1_dt_c = zeros(1, n);
    dQ2_dt_c = zeros(1, n);
    dQ3_dt_c = zeros(1, n);
    Q1_c = zeros(1, n);
    Q2_c = zeros(1, n);
    Q3_c = zeros(1, n);

    for i = 2:n
        dQ1_dt_c(i) = -(F1_p(i) - F1_p(i-1)) / dx;
        dQ2_dt_c(i) = -(F2_p(i) - F2_p(i-1)) / dx + (1 / gamma) * (rho(i) * T(i)) * ((A(i) - A(i-1)) / dx);
        dQ3_dt_c(i) = -(F3_p(i) - F3_p(i-1)) / dx;

        Q1_c(i) = Q1_p(i) + dt * dQ1_dt_c(i);
        Q2_c(i) = Q2_p(i) + dt * dQ2_dt_c(i);
        Q3_c(i) = Q3_p(i) + dt * dQ3_dt_c(i);
    end
end
