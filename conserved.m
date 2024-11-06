function [Q1, Q2, Q3] = conserved(rho, A, v, T, gamma)
    Q1 = rho.*A;
    Q2 = rho.*A.*v;
    Q3 = rho .* A .* (T / (gamma - 1) + 0.5*gamma* v.^2);
end

