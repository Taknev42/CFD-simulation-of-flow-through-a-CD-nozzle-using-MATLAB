function [rho, v, T, P, Mach_no] = update_primitives(Q1, Q2, Q3, A, gamma)
    rho = Q1 ./ A;
    v = Q2 ./ Q1;
    T = (gamma - 1) * (Q3 ./ Q1 - 0.5 * gamma * v.^2);
    P = rho .* T;
    Mach_no = v ./ sqrt(T);
end