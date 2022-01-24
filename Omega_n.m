function [price] = Omega_n(u_n,T,eta,kappa_n)
%Returns the price of electricity given the parameters and power use

p_n = u_n(1:T);
price = eta'*p_n + kappa_n*(p_n'*p_n);

end

