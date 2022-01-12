function [price] = Omega_n(p_n,eta,kappa_n)
%Returns the price of electricity given the parameters and power use

price = eta'*p_n + kappa_n*(p_n'*p_n);

end

