%{
Testing the IEEE 13 node feeder.
Jan 2022 - Abhishek B.
%}
set(0,'DefaultFigureWindowStyle','docked')
%% REES model parameters and such

% Customers, time steps, supply point numbers
N = 3; 
T = 3;
K = 3;


arrival = 0; 
departure = T; 

capacities = [42,21,24,20,30,16,23,16.5,28,60,90,75];
b = zeros(1,N); 
for i=1:N
    b(i) = capacities(randi([1,12]));
end
s_bar = 10; 
sig_hat = (0.3+0.2*rand(1,N)).*b;
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
sig_star = zeros(1,N); 
for i=1:N
    sig_star(i) = ceil(sig_hat(i) + (sig_u(i)-sig_hat(i))*rand);
end
p_u = 6.6;  % max charge rate
p_l = -6.6; % min charge rate

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

eta_simple = eta(1:16:end);

% Regularization paramters for cost
kappa = 5e-4; %zeros(1,N);

%% Unbalanced dist grid model

% Matrix Theta - this has size 29x600

% This has to be manually made each time anyway, so no reason to try and
% automate it at this stage.
% theta = zeros(K,N);
theta = eye(3,3);

% R X matrices
% miles to feet
mi2ft = 1/5280;
% Config 601
R_601 = [0.3465, 0.1560, 0.1580;
     0.1560, 0.3375, 0.1535;
     0.1580, 0.1535, 0.3413]; 
X_601 = [1.0179, 0.5017, 0.4236;
     0.5017, 1.0478, 0.3849;
     0.4236, 0.3849, 1.0348];
Z_601 = (R_601+1i*X_601)*mi2ft;
omega = exp(-2*pi*1i/3);
for i=1:3
    for j=1:3        
        R(i,j) = 2*real(conj(Z_601(i,j))*(omega^(i-j)));   %2*real(Z_601(1,1)*2000);
        X(i,j) = -2*imag(conj(Z_601(i,j))*(omega^(i-j)));  %2*imag(Z_601(1,1)*2000);
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


%% Old remnant code for baseline voltages
 

v_temp = [1.05, 0.95, 1]; 
V_baseline = zeros(K,T);
for j=1:T
    for i=1:K
        V_baseline(i,j) = [v_temp(j)-0.01*v_temp(j)] + 0.02*rand;
    end
end

% figure()
% for i=1:K
%     plot(V_baseline(i,:), 'k-o')
%     ylim([0.9,1.1])
%     hold on
%     grid on
%     set(gca,'ytick',[0.0:0.01:1.1])
%     plot(1.046*ones(1,T), 'r-')
%     plot(0.954*ones(1,T), 'r-')
% end
% disp(V_baseline)

v_upper = 1.046*ones(K*T,1);
v_lower = 0.954*ones(K*T,1);

w = [v_upper - V_baseline(:);
     -v_lower + V_baseline(:)];
%% Baseline voltages and vec w
% 
% load('baseload.mat');
% 
% v_upper = 1.046*ones(K*T,1);
% v_lower = 0.954*ones(K*T,1);
% 
% V_base_stacked = V_base(:);
% 
% w = [v_upper - V_base_stacked;
%      -v_lower - V_base_stacked];

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

% %% Costs
% 
% costs={};
% for i=1:N
%     costs{i} = @(u) eta_simple'*u(1:T) + kappa*u(1:T)'*u(1:T); %Omega_n(u, T, eta, kappa); % + Indicator_n(u, T, w, s_bar(i), mu_u, mu_l, a(i), d(i), p_u, p_l, alpha_u(i), alpha_l(i), e(i));
% end


%% ADMM

iteration_limit=1000;

step_size = 0.0001;
err_tol = 1e-6;

% initialise storage
lambda_step = {};
nu_step = {};
u_step = {};

% initialise the variables
u = zeros(2*T+2*K*T,N);
lambda = zeros(6*T,N);
nu = zeros(6*T,N);

iter = 0;
err = 1;
err_step = [err];

while (err>err_tol) && (iter<iteration_limit)
    % Update u and lambda
    iter = iter+1;
    disp(iter)
    
    u_old = u;
    lambda_old = lambda;
    nu_old = nu;
    
    
    parfor i = 1:N
        % initialise sdpvar
        
        x = sdpvar(2*T+2*K*T,1);
        p = x(1:T);
        q = x(T+1:2*T);
        s = x(2*T+1:end);
        
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
        constraints = [constraints, ones(1,T)*p == e(i)];
        
        % slack variable constraint
        constraints = [constraints, zeros(2*K*T,1) <= s];
        
        % objective
        
        % Building the sum of self state and neighbors
        % states
        temp = 0;
        for j=1:length(neighbors{i})
            temp = temp +( lambda_old(:,i) + lambda_old(:,neighbors{i}(j)) );
        end
        
        obj = (eta_simple'*p + kappa*(p'*p)) +...
              (step_size/(n_deg(i)*4))*norm( (1/step_size)*(Psi{i}*x - (1/N)*w) - (1/step_size)*nu_old(:,i) +  temp)^2;
        
        options = sdpsettings('verbose',0);
        optimize(constraints, obj, options);
        u(:,i) = value(x); 
        
        lambda(:,i) = (1/(2*n_deg(i)))*( temp - (1/step_size)*nu_old(:,i) + (1/step_size)*( Psi{i}*u(:,i) - (1/N)*w )  ); 
    end 
    
    u_step{iter} = u;
    lambda_step{iter} = lambda;
    
    % Updating nu
    parfor i=1:N
        temp1 = 0;
        for j=1:length(neighbors{i})
            temp1 = temp1 + lambda(:,i) - lambda(:,j);
        end
        
        nu(:,i) = nu_old(:,i) + step_size*temp1;
    end
    
    u_err = norm(u_old-u,'fro')^2;
    nu_err = norm(nu_old-nu,'fro')^2;
    lambda_err = norm(lambda_old - lambda, 'fro')^2;
    
    err=max([u_err, nu_err, lambda_err]);
    err_step = [err_step, err];
end

figure()
plot(err_step(3:end))

%% Reconstruct Voltages from solutions

% disp(V_baseline)

V_soln = zeros(K*T,1);

temp = zeros(K*T,1);
for i=1:N
    temp = temp + D_stack{i}*u(1:T, i) + E_stack{i}*u(T+1:2*T, i);
end

V_soln = V_baseline(:) + temp;

V_soln_reshaped = reshape(V_soln, [K,T]);

%% Plots
figure()
for i=1:K
    plot(V_baseline(i,:), 'rx--')
    ylim([0.9,1.1])
    hold on
    plot(1.046*ones(1,T), 'k-')
    hold on
    plot(0.954*ones(1,T), 'k-')
    grid on
    set(gca,'ytick',[0.0:0.01:1.1])
    hold on
    plot(V_soln_reshaped(i,:),'go-')
end

