%{
Testing the IEEE 13 node feeder.
Jan 2022 - Abhishek B.
%}
set(0,'DefaultFigureWindowStyle','docked')
rng('default')
seed = rng;
%% REES model parameters and such

% Customers, time steps, supply point numbers
cust_at_supply = 25;
N = 3*cust_at_supply; 
T = 48;
K = 3;


% Customer Parameters
arrival = 0; 
departure = T; 

capacities = [42,21,24,20,30,16,23,16.5,28,60,90,75]; %batteries in kWh
b = zeros(1,N); 
for i=1:N
    b(i) = capacities(randi([1,12]));
end
s_bar = 10; % Inverter capacity in  
sig_hat = (0.3+0.2*rand(1,N)).*b;
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
sig_star = zeros(1,N); 
for i=1:N
    sig_star(i) = ceil(sig_hat(i) + (sig_u(i)-sig_hat(i))*rand);
end
p_u = 6.6;  % max charge rate in kWATTS
p_l = -6.6; % min charge rate in kWATTS

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

theta = [ones(1,cust_at_supply), zeros(1,cust_at_supply), zeros(1,cust_at_supply);
         zeros(1,cust_at_supply), ones(1,cust_at_supply), zeros(1,cust_at_supply);
         zeros(1,cust_at_supply), zeros(1,cust_at_supply), ones(1,cust_at_supply)];
     
% R X matrices

mi2ft = 1/5280;  % miles to feet
R_601 = [0.3465, 0.1560, 0.1580;
         0.1560, 0.3375, 0.1535;
         0.1580, 0.1535, 0.3413];
X_601 = [1.0179, 0.5017, 0.4236;
         0.5017, 1.0478, 0.3849;
         0.4236, 0.3849, 1.0348];
Z_601 = (R_601+1i*X_601)*mi2ft*2000;

R = zeros(K,K);
X = zeros(K,K);
omega = exp(-2*pi*1i/3);

for i=1:K
    for j=1:K
        R(i,j) = 2*real( conj(Z_601(i,j))*(omega^(i-j))  );
        X(i,j) = -2*imag( conj(Z_601(i,j))*(omega^(i-j))  );
    end
end


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

% Baseline Voltage in kV.
v0=2.4;

v_upper = ((1.05)^2)*ones(K*T,1);
v_lower = ((0.95)^2)*ones(K*T,1);

U = U_base(1:3, :);
U_stacked = U(:);

V_b = V_base(1:3,:).^2;
V_stacked = V_b(:);

w = [v_upper - V_stacked;
     -v_lower - V_stacked];

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

%% ADMM

iteration_limit=100;

step_size = 0.0001;
progress_tol = 1e-3;

% initialise storage
lambda_ADMM_step = {};
nu_ADMM_step = {};
u_ADMM_step = {};

% initialise the variables
u_ADMM = zeros(2*T+2*K*T,N);
lambda_ADMM = zeros(2*K*T,N);
nu_ADMM = zeros(2*K*T,N);

iter = 0;
change = 1;
U_ADMM_err = [];
Lambda_ADMM_err = [];
Nu_ADMM_err = [];
Change_ADMM_step = [change];

while ((change>progress_tol) || isnan(change)) && (iter<iteration_limit) 
    % Update u and lambda
    iter = iter+1;
    disp(strcat("iteration - ",int2str(iter)," | progress - ", num2str(change)))
    
    u_ADMM_old = u_ADMM;
    lambda_ADMM_old = lambda_ADMM;
    nu_ADMM_old = nu_ADMM;
    
    
    parfor i = 1:N
        % initialise sdpvar
        
%         disp('working 1')
        x = sdpvar(2*T+2*K*T,1);
        p = x(1:T);
        q = x(T+1:2*T);
        s = x(2*T+1:end);
        
%         disp('working 2')
        % constraints
        constraints = [];
        
        % inverter capacity constraint
        constraints = [constraints, p.^2 + q.^2 <= (s_bar^2)*ones(T,1)];
        
        % over/undercharge prevention constraints
        constraints = [constraints, tril(ones(T,T))*p <= alpha_u(i)];
        constraints = [constraints, alpha_l(i) <= tril(ones(T,T))*p];
        
        % charge discharge limit constraints
        constraints = [constraints, p <= p_u*ones(T,1)];
        constraints = [constraints, p_l*ones(T,1) <= p];
        
        % charge demand constraint
        constraints = [constraints, ones(1,T)*p >= e(i)];
        
        % slack variable constraint
        constraints = [constraints, zeros(2*K*T,1) <= s];
        
        % objective
        
        % Building the sum of self state and neighbors
        % states
        temp = 0;
        for j=1:length(neighbors{i})
            temp = temp +( lambda_ADMM_old(:,i) + lambda_ADMM_old(:,neighbors{i}(j)) );
        end
        
        obj = (eta'*p + kappa*(p'*p)) +...
              (step_size/(n_deg(i)*4))*norm( (1/step_size)*(Psi{i}*(x/1000) - (1/N)*w) - (1/step_size)*nu_ADMM_old(:,i) +  temp)^2;
        
        options = sdpsettings('verbose',0);
        optimize(constraints, obj, options);
        u_ADMM(:,i) = value(x)/1000; 
        
        lambda_ADMM(:,i) = (1/(2*n_deg(i)))*( temp - (1/step_size)*nu_ADMM_old(:,i) + (1/step_size)*( Psi{i}*u_ADMM(:,i) - (1/N)*w )  ); 
    end 
    
    u_ADMM_step{iter} = u_ADMM;
    lambda_ADMM_step{iter} = lambda_ADMM;
    
    % Updating nu
    parfor i=1:N
        temp1 = 0;
        for j=1:length(neighbors{i})
            temp1 = temp1 + lambda_ADMM(:,i) - lambda_ADMM(:,j);
        end
        
        nu_ADMM(:,i) = nu_ADMM_old(:,i) + step_size*temp1;
    end
    
    nu_ADMM_step{iter} = nu_ADMM;
    
    u_ADMM_err = norm(u_ADMM_old-u_ADMM,'fro')^2;
    nu_ADMM_err = norm(nu_ADMM_old-nu_ADMM,'fro')^2;
    lambda_ADMM_err = norm(lambda_ADMM_old - lambda_ADMM, 'fro')^2;
    
    U_ADMM_err = [U_ADMM_err, u_ADMM_err];
    Lambda_ADMM_err = [Lambda_ADMM_err, lambda_ADMM_err];
    Nu_ADMM_err = [Nu_ADMM_err, nu_ADMM_err];
    
    change=max([u_ADMM_err, nu_ADMM_err, lambda_ADMM_err]);
    Change_ADMM_step = [Change_ADMM_step, change];
end


%% Censored ADMM Solution


iteration_limit=100;

step_size = 0.0001;
progress_tol = 1e-3;

% Censoring parameters
censoring_rho = 0.7;
censoring_alpha = 0.5;

% initialise storage
lambda_cens_step = {};
nu_cens_step = {};
u_cens_step = {};

% initialise the variables
u_cens = zeros(2*T+2*K*T,N);
lambda_cens = zeros(2*K*T,N);
nu_cens = zeros(2*K*T,N);

iter = 0;
change_cens = 1;
U_cens_err = [];
Lambda_cens_err = [];
Nu_cens_err = [];
Change_cens_step = [change_cens];

while ((change_cens>progress_tol) || isnan(change_cens)) && (iter<iteration_limit) 
    % Update u and lambda
    iter = iter+1;
    disp(strcat("iteration - ",int2str(iter)," | progress - ", num2str(change_cens)))
    
    u_cens_old = u_cens;
    lambda_cens_old = lambda_cens;
    nu_cens_old = nu_cens;
    
    
    parfor i = 1:N
        % initialise sdpvar
        
%         disp('working 1')
        x = sdpvar(2*T+2*K*T,1);
        p = x(1:T);
        q = x(T+1:2*T);
        s = x(2*T+1:end);
        
%         disp('working 2')
        % constraints
        constraints = [];
        
        % inverter capacity constraint
        constraints = [constraints, p.^2 + q.^2 <= (s_bar^2)*ones(T,1)];
        
        % over/undercharge prevention constraints
        constraints = [constraints, tril(ones(T,T))*p <= alpha_u(i)];
        constraints = [constraints, alpha_l(i) <= tril(ones(T,T))*p];
        
        % charge discharge limit constraints
        constraints = [constraints, p <= p_u*ones(T,1)];
        constraints = [constraints, p_l*ones(T,1) <= p];
        
        % charge demand constraint
        constraints = [constraints, ones(1,T)*p >= e(i)];
        
        % slack variable constraint
        constraints = [constraints, zeros(2*K*T,1) <= s];
        
        % objective
        
        % Building the sum of self state and neighbors
        % states
        temp = 0;
        for j=1:length(neighbors{i})
            temp = temp +( lambda_cens_old(:,i) + lambda_cens_old(:,neighbors{i}(j)) );
        end
        
        obj = (eta'*p + kappa*(p'*p)) +...
              (step_size/(n_deg(i)*4))*norm( (1/step_size)*(Psi{i}*(x/1000) - (1/N)*w) - (1/step_size)*nu_cens_old(:,i) +  temp)^2;
        
        options = sdpsettings('verbose',0);
        optimize(constraints, obj, options);
        u_ADMM(:,i) = value(x)/1000; 
        
        lambda_ADMM(:,i) = (1/(2*n_deg(i)))*( temp - (1/step_size)*nu_cens_old(:,i) + (1/step_size)*( Psi{i}*u_cens(:,i) - (1/N)*w )  ); 
    end 
    
    u_cens_step{iter} = u_cens;
    lambda_cens_step{iter} = lambda_cens;
    
    % Updating nu
    parfor i=1:N
        temp1 = 0;
        for j=1:length(neighbors{i})
            temp1 = temp1 + lambda_cens(:,i) - lambda_cens(:,j);
        end
        
        nu_cens(:,i) = nu_cens_old(:,i) + step_size*temp1;
    end
    
    nu_cens_step{iter} = nu_cens;
    
    u_cens_err = norm(u_cens_old-u_cens,'fro')^2;
    nu_cens_err = norm(nu_cens_old-nu_cens,'fro')^2;
    lambda_cens_err = norm(lambda_cens_old - lambda_cens, 'fro')^2;
    
    U_cens_err = [U_cens_err, u_cens_err];
    Lambda_cens_err = [Lambda_cens_err, lambda_cens_err];
    Nu_cens_err = [Nu_cens_err, nu_cens_err];
    
    change_cens=max([u_cens_err, nu_cens_err, lambda_cens_err]);
    Change_cens_step = [Change_cens_step, change_cens];
end



%% Reconstruct Voltages from solutions

% disp(V_baseline)

V_soln = zeros(K*T,1);

temp = zeros(K*T,1);
for i=1:N
    temp = temp + D_stack{i}*u_ADMM(1:T, i) + E_stack{i}*u_ADMM(T+1:2*T, i);
end

% temp=temp/1000;

V_soln = V_stacked + (temp/(v0^2));

V_soln_reshaped = reshape(V_soln, [K,T]);

%% Plots
figure()
for i=1:K
    plot(sqrt(V_b(i,:)), 'rx--')
    hold on
    plot(1.05*ones(1,T), 'k-')
    hold on
    plot(0.95*ones(1,T), 'k-')
    grid on
    hold on
    plot(sqrt(V_soln_reshaped(i,:)),'go-')
end


