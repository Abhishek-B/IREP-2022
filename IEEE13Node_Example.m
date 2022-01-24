%{
Testing the IEEE 13 node feeder.
Jan 2022 - Abhishek B.
%}
set(0,'DefaultFigureWindowStyle','docked')
%% REES model parameters and such

% Customers, time steps, supply point numbers
N = 600;
T = 48;
K = 29;

% power profile vectors (empty for now)
p = zeros(T,N);
q = zeros(T,N);
s = zeros(2*K*T, N);
u = [p;q;s];

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

% Matrix Theta - this has size 29x600

theta = zeros(29, 600);

% The only way I can make sense of the numbers in Nanduni's paper is if
% each supply point has its own number of customers. So each customer is
% essentially only connected to a single phase...
% yeah this incrementation and labelling of customers is a bit confusing...

theta(4,1:26) = 1;
theta(5,27:52) = 1;
theta(6,53:78) = 1;
theta(7,79:104) = 1;
theta(8,105:130) = 1;
theta(9,131:156) = 1;
theta(10,157:182) = 1;
theta(11, 183:208) = 1;
theta(12, 209:234) = 1;
theta(13, 235:260) = 1;
theta(17, 261:286) = 1;
theta(18, 287:312) = 1;
theta(19, 313:338) = 1;
theta(20, 339:364) = 1;
theta(21, 365:390) = 1;
theta(22, 391:416) = 1;
theta(23, 417:442) = 1;
theta(24, 443:468) = 1;
theta(25, 469:495) = 1; % Supply point 10 which has 27 customers
theta(26, 496:522) = 1; % Supply point 11 which has 27 customers
theta(27, 523:548) = 1;
theta(28, 549:574) = 1;
theta(29, 575:600) = 1;

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

Psi = {};
for i=1:N
    Psi{i} = [Gamma{i}, Xi{i}, eye(2*K*T)];
end

%% Baseline voltages and vec w

load('baseload.mat');

v_upper = 1.046*ones(K*T,1);
v_lower = 0.954*ones(K*T,1);

w = [v_upper - V_base;
     -v_lower - V_base];

%% Comms network

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

%% Costs

costs={};
for i=1:N
    costs{i} = @(u) Omega_n(u, T, eta, kappa); % + Indicator_n(u, T, w, s_bar(i), mu_u, mu_l, a(i), d(i), p_u, p_l, alpha_u(i), alpha_l(i), e(i));
end

%% ADMM

iteration = 1000;
lambda = zeros(2*K*T,N);       % primal value for each customer
lambda_state = zeros(2*K*T,N); % state value for each customer, the neighbors states
                      % can be found from different columns of this.

lambda_step = {}; % cell to hold each iteration value of the primal
nu_step = {}; % cell to hold each iteration value of the dual

xi = zeros(2*K*T,N);      % the difference between primal value and state value
nu = zeros(2*K*T,N);  % dual variabels for each customer

transmission = zeros(N,iteration);

for k=1:iteration
    
    % Update of the primal variable
    for i=1:N
        % Building the sum of difference between self state and neighbors
        % states
        temp=0;
        for j=1:length(neighbors{i})
            temp = temp +( lambda_state(:,i) + lambda_state(:,neighbors{i}(j)) );
        end
        
        
        
        
        % the overall cost function for the primal update
        fun = @(x)costs{i}(x) + x'*(lambda(:,i) - c*temp) + c*n_deg(i)*norm(x)^2;
        
        % updating the primal through a subminimization
        lambda(:,i) = fminsearch(fun, zeros(p,1));
        
        
        
        
        
        
        % Computing the difference between current state and primal update
        xi(:,i) = lambda_state(:,i) - lambda(:,i); 
        %disp(xi(:,i)')
        % Transmission loop
        if ( norm(xi(:,i))^2 - alpha*(rho^k) )>=0
            %disp(norm(xi(:,i))^2 - alpha*(rho^k))
            %disp("true")
            lambda_state(:,i) = lambda(:,i);
            transmission(i,k) = 1;
        else
            %disp("false")
        end
    end
    
    % Dual variable update
    for i=1:n
        
        % building the sum of the differences between self state and
        % neighbors current states
        temp2 = 0;
        for j=1:length(neighbors{i})
            temp2 = temp2 + ( lambda_state(:,i) - lambda_state(:,neighbors{i}(j)) );
        end
        
        % Updating the dual variable
        nu(:,i) = nu(:,i) + c*temp2;
    end
    
    % Storing the primal and dual updates for the iteration
    lambda_step{k} = lambda;
    nu_step{k} = nu;
end











 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
% % baseline voltage
% % I'm gonna set this up myself by choosing random values. The dataset from
% % ausgrid doesn't seem to have the link up anymore, and i wouldn't know how
% % to convert nodal voltage data to the p.u. system. 
% 
% 
% % these have 49 values and not 48...so maybe the time steps are 49=T
% x = [   1,    4,     5, 8,   11,    16,   17,    20,   23,   32,    34,    37,   40, 42,   45, 48];
% xq = [1:1:48];
% v = [1.04, 1.01, 1.023, 1, 0.95, 0.955, 0.96, 0.974, 0.98, 0.98, 0.976, 0.974, 0.97,  1, 1.01,  1];
% 
% v_temp = interp1(x,v,xq);
% 
% % figure()
% % plot(xq, v_temp, 'r-o')
% % ylim([0.9,1.1])
% 
% V_baseline = zeros(K,T);
% 
% for j=1:T
%     for i=1:K
%         V_baseline(i,j) = [v_temp(j)-0.01*v_temp(j)] + 0.02*rand;
%     end
% end
% 
% % figure()
% % for i=1:29
% %     plot(xq, V_baseline(i,:), 'k-o')
% %     ylim([0.9,1.1])
% %     hold on
% %     grid on
% %     set(gca,'ytick',[0.0:0.01:1.1])
% % end
% % plot(xq, 1.046*ones(1,48), 'r-')
% % plot(xq, 0.954*ones(1,48), 'r-')
% 
% % The plot seem okay-ish for these baseline voltages. Eyeballing it, they
% % seem roughly similar to Nanduni's plot, so I'll go with this.
% 
% V_base = [];
% for i=1:T
%     V_base = [V_base; V_baseline(:,i)];
% end
% 
% v_upper = 1.046*ones(K*T,1);
% v_lower = 0.954*ones(K*T,1);
% 
% w = [v_upper - V_base;
%      -v_lower - V_base];

