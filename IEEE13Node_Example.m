%{
Testing the IEEE 13 node feeder.
Jan 2022 - Abhishek B.
%}

%% REES model parameters and such

% Customers and time steps
N = 600;
T = 48;

% power profile vectors (empty for now)
p = zeros(T,N);
q = zeros(T,N);

% Stacking design parameters into row vectors, each elements corresponds to
% customer i
% setting to zero for now, will fill in soon

a = zeros(1,N); % arrival time
d = zeros(1,N); % departure time
b = zeros(1,N); % battery capacity - this isn't used in the manuscript anywhere
s_bar = zeros(1,N); % inverter capacity
sig_hat = (0.3+0.2*rand(1,N)).*b; % initial SoC
sig_star = zeros(1,N); % target SoC
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
p_u = 6.6;  % max charge rate
p_l = -6.6; % min charge rate
mu_u = 0.9; % charge efficiency
mu_l = 1.1; % discharge efficiency


% variables defined for linear systems
alpha_l = sig_l - sig_hat;
alpha_u = sig_u - sig_hat;
e = sig_star - sig_hat;

%% Operational cost parameters

% price column vector
eta = zeros(T,1);

% Regularization paramters for cost
kappa = 5e-4; %zeros(1,N);

%% Unbalanced dist grid model


