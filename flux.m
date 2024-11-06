function [F1, F2, F3] = flux(Q1, Q2, Q3, gamma)
    F1 = Q2;
    F2 = (Q2.^2)./(Q1) + ((gamma-1)/gamma)*(Q3 - 0.5*gamma*(Q2.^2)./(Q1));
    F3 = gamma*(Q2.*Q3./Q1) - 0.5*(gamma*(gamma-1))*((Q2.^3)./(Q1.^2)); 
end
