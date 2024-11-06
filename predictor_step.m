function [Q1_p, Q2_p, Q3_p, dQ1_dt_p, dQ2_dt_p, dQ3_dt_p] = predictor_step(Q1, Q2, Q3, F1, F2, F3, dx, dt, A, rho, T, gamma)
    n = length(Q1);
    dQ1_dt_p = zeros(1, n);
    dQ2_dt_p = zeros(1, n);
    dQ3_dt_p = zeros(1, n);
    Q1_p = zeros(1, n);
    Q2_p = zeros(1, n);
    Q3_p = zeros(1, n);

    for i = 1:n-1
        dQ1_dt_p(i) = -(F1(i+1) - F1(i)) / dx;
        dQ2_dt_p(i) = -(F2(i+1) - F2(i)) / dx + (1 / gamma) * (rho(i) * T(i)) * ((A(i+1) - A(i)) / dx);
        dQ3_dt_p(i) = -(F3(i+1) - F3(i)) / dx;

        Q1_p(i) = Q1(i) + dt * dQ1_dt_p(i);
        Q2_p(i) = Q2(i) + dt * dQ2_dt_p(i);
        Q3_p(i) = Q3(i) + dt * dQ3_dt_p(i);
    end
end