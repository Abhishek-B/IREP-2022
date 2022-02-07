%{
Simplified model for the IREP stuff.
-Abhishek B.
%}

set(0,'DefaultFigureWindowStyle','docked')
%% REES model parameters and such

% Customers, time steps, supply point numbers
N = 23; 
T = 48;
K = 29;

% Changing these to simplest possible scenario.
arrival = 0; 
departure = T; 

capacities = [42,21,24,20,30,16,23,16.5,28,60,90,75];
b = zeros(1,N); 
for i=1:N
    b(i) = capacities(randi([1,12]));
end
s_bar = 10; % inverter capacity
sig_hat = (0.3+0.2*rand(1,N)).*b; % initial SoC
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
sig_star = zeros(1,N); % target SoC - this isnt listed in the paper, so I 
                       % will set this to between sig_hat and sig_upper,
                       % rounded to the nearest integer
for i=1:N
    sig_star(i) = ceil(sig_hat(i) + (sig_u(i)-sig_hat(i))*rand);
end
p_u =  6.6; % max charge rate
p_l = -6.6; % min charge rate
q_u =  6.6; % Simple box constraints for p and q, instead of nonlinear
q_l = -6.6;

% variables defined for linear systems
alpha_l = sig_l - sig_hat;
alpha_u = sig_u - sig_hat;
e = sig_star - sig_hat;

%% Operational cost parameters

% price column vector
eta = zeros(T,1);
eta(1:4) = 0.23;
eta(5:16)=0.54;
eta(17:20) = 0.23;
eta(21:38) = 0.15;
eta(39:48) = 0.23;

% Regularization paramters for cost
kappa = 5e-4; %zeros(1,N);

%% Unbalanced dist grid model

% Matrix Theta - this has size 29x600

theta = zeros(29, N);
theta(4, 1) = 1;
theta(5, 2) = 1;
theta(6, 3) = 1;
theta(7, 4) = 1;
theta(8, 5) = 1;
theta(9, 6) = 1;
theta(10, 7) = 1;
theta(11, 8) = 1;
theta(12, 9) = 1;
theta(13, 10) = 1;
theta(17, 11) = 1;
theta(18, 12) = 1;
theta(19, 13) = 1;
theta(20, 14) = 1;
theta(21, 15) = 1;
theta(22, 16) = 1;
theta(23, 17) = 1;
theta(24, 18) = 1;
theta(25, 19) = 1; 
theta(26, 20) = 1; 
theta(27, 21) = 1;
theta(28, 22) = 1;
theta(29, 23) = 1;

% Loading R and X
load('RX.mat');

% matrices D and E
D = -R*theta;
E = -X*theta;

D_stack = {};
for i=1:N
    temp = D(:,i);
    for j=1:(T-1)
        temp = blkdiag(temp,D(:,i));
    end
    D_stack{i} = temp;
end

E_stack = {};
for i=1:N
    temp = E(:,i);
    for j=1:(T-1)
        temp = blkdiag(temp,E(:,i));
    end
    E_stack{i} = temp;
end


% Matrices Gamma and Xi
Gamma = {};
Xi = {};
for i=1:N
    Gamma{i} = [ D_stack{i} ;
                -D_stack{i} ];
    Xi{i} = [ E_stack{i} ;
             -E_stack{i} ];
end

%% matrices for ADMM

A = {};
B = {};
M = tril(ones(T,T));
One = ones(T,1);

for i=1:N
    Id_stack = [];
    for j=1:N
        if j==i
            Id_stack = [Id_stack; eye(T,T)];
        else
            Id_stack = [Id_stack; zeros(T,T)];
        end
    end
    M_stack = [];
    for j=1:N
        if j==i
            M_stack = [M_stack; M];
        else
            M_stack = [M_stack; zeros(T,T)];
        end
    end
    One_stack = [];
    for j=1:N
        if j==i
            One_stack = [One_stack; One'];
        else
            One_stack = [One_stack; zeros(1,T)];
        end
    end
    zero_stack = zeros(N*T, T);
    
    A{i} = [  Gamma{i}  ;
              Id_stack  ;
             -Id_stack  ;
              M_stack   ;
             -M_stack   ;
             -One_stack ;
             zero_stack ;
             zero_stack ];
    
    B{i} = [    Xi{i}   ;
             zero_stack ;
             zero_stack ;
             zero_stack ;
             zero_stack ;
             zeros(N,T) ;
              Id_stack  ;
             -Id_stack  ];
end




Psi = {};
nrows = size(A{1},1);
for i=1:N
    Psi{i} = [A{i}, B{i}, eye(nrows,nrows)];
end

%% Old remnant code for baseline voltages
 
x = [   1,    4,     5, 8,   11,    16,   17,    20,   23,   32,    34,    37,   40, 42,   45, 48];
xq = [1:1:48];
v = [1.04, 1.01, 1.023, 1, 0.95, 0.955, 0.96, 0.974, 0.98, 0.98, 0.976, 0.974, 0.97,  1, 1.01,  1];

v_temp = interp1(x,v,xq);

V_baseline = zeros(K,T);

v_upper = 1.046*ones(K*T,1);
v_lower = 0.954*ones(K*T,1);

w = [v_upper - V_baseline(:);
     -v_lower - V_baseline(:)];
 
%% Updating w into W for admm matrices

W = [w];
% p_upper and p_lower are the same for everyone, so this next part combines
% all the p constraints


W = [w; p_u*ones(N*T,1); -p_l*ones(N*T,1)];

% adding the alpha constraints

for i=1:N
    W = [ W ; alpha_u(i)*ones(T,1) ];
end

for i=1:N
    W = [ W ; -alpha_l(i)*ones(T,1) ];
end

% adding e(i) constraint

for i=1:N
    W = [ W ; -e(i)];
end

% adding q constraint - q_upper and lower are the same so this is like the
% p thing

W = [W ; q_u*ones(N*T,1) ; -q_l*ones(N*T,1)];

%% Graph parameters setup
% Comms network
Adjacency = ones(N,N) - eye(N,N);
G = graph(Adjacency);
edge = G.Edges;
n_deg = degree(G);
neighbors = {};
for i=1:N
    neighbors{i} = find(Adjacency(:,i));
end


