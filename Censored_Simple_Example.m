%{
Testing simple censored example.
Feb 2022 - Abhishek B.
%}
set(0,'DefaultFigureWindowStyle','docked')
rng('default')
seed = rng;
%% REES model parameters and such

% Customers, time steps, supply point numbers
cust_at_supply = 50;
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
s_bar = 12; % Inverter capacity in  
sig_hat = (0.3+0.2*rand(1,N)).*b;
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
sig_star = zeros(1,N); 
for i=1:N
    sig_star(i) = ceil(sig_hat(i) + (sig_u(i)-sig_hat(i))*rand);
end
p_u = 7;  % max charge rate in kWATTS
p_l = -7; % min charge rate in kWATTS

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
D = -(1/1000)*R*theta;
E = -(1/1000)*X*theta;

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

percentage = 0.05;

v_upper = ((1+percentage)^2)*ones(K*T,1);
v_lower = ((1-percentage)^2)*ones(K*T,1);

U = U_base(1:3, :);
U_stacked = U(:);

V_b = V_base(1:3,:).^2;
V_stacked = V_b(:);

w = [v_upper - V_stacked;
     -v_lower + V_stacked];

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

%% ADMM Parameters
iteration_limit=30;

step_size = 100;
progress_tol = 1e-3;
%% ADMM

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
time_ADMM = [];

disp("Computing ADMM solutions")

while ((change>progress_tol) || isnan(change)) && (iter<iteration_limit) 
    % Update u and lambda
    iter = iter+1;
    tic;
    
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
              (step_size/(n_deg(i)*4))*norm( (1/step_size)*(Psi{i}*(x) - (1/N)*w) - (1/step_size)*nu_ADMM_old(:,i) +  temp)^2;
        
        options = sdpsettings('verbose',0);
        optimize(constraints, obj, options);
        u_ADMM(:,i) = value(x); 
        
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

    time_taken = toc;
    time_ADMM = [time_ADMM, time_taken];

    disp(strcat("iter - ",int2str(iter)," | prog - ", num2str(round(change,11)), " | time - ", num2str(round(time_taken/60,1))))

end

%% ADMM Plot
% 
% % disp(V_baseline)
% 
% V_soln = zeros(K*T,1);
% 
% p_ADMM = {};
% q_ADMM = {};
% % p_cens = {};
% % q_cens = {};
% for i=1:N
%     p_ADMM{i} = u_ADMM(   1:T   , i);
%     q_ADMM{i} = u_ADMM( T+1:2*T , i);
% %     p_cens{i} = u_cens(   1:T   , i);
% %     q_cens{i} = u_cens( T+1:2*T , i);
% end
% 
% V_control_ADMM = zeros(K*T,1);
% % V_control_cens = zeros(K*T,1);
% 
% for i=1:N
%     V_control_ADMM = V_control_ADMM + D_stack{i}*p_ADMM{i} + E_stack{i}*q_ADMM{i};
% %     V_control_cens = V_control_cens + D_stack{i}*p_cens{i} + E_stack{i}*q_cens{i};
% end
% 
% V_soln_ADMM = V_stacked + (V_control_ADMM/(v0^2));
% % V_soln_cens = V_stacked + (V_control_cens/(v0^2));
% 
% V_ADMM_reshaped = reshape(V_soln_ADMM, [K,T]);
% % V_cens_reshaped = reshape(V_soln_cens, [K,T]);
% 
% p_aggregate = {};
% q_aggregate = {};
% 
% for i=1:K
%     p = zeros(48,1);
%     q = zeros(48,1);
%     if i==1
%         for j=1:30
%             p = p + p_ADMM{j};
%             q = q + q_ADMM{j};
%         end
%     elseif i==2
%         for j=31:60
%             p = p + p_ADMM{j};
%             q = q + q_ADMM{j};
%         end
%     elseif i==3
%         for j=61:90
%             p = p + p_ADMM{j};
%             q = q + q_ADMM{j};
%         end
%     end
%     p_aggregate{i} = p;
%     q_aggregate{i} = q;
% end
% 
% for i=1:K
%     figure()
%     tiledlayout(3,1)
%     
%     ax1 = nexttile;
%     plot(ax1,sqrt(V_b(i,:)), 'rx--')
%     hold on
%     plot(ax1,(1+percentage)*ones(1,T), 'k-')
%     hold on
%     plot(ax1,(1-percentage)*ones(1,T), 'k-')
%     hold on
%     plot(ax1,sqrt(V_ADMM_reshaped(i,:)),'g--')
%     hold on
% %     plot(ax1,sqrt(V_cens_reshaped(i,:)),'bo')
%     grid on
%     title(ax1, strcat("Baseline Voltage and Control Voltage of Node ", num2str(i)))
%     
%     ax2 = nexttile;
%     plot(ax2,p_aggregate{i})
%     hold on
%     plot(ax2,zeros(48,1),'k--')
%     grid on
%     title(ax2, strcat("Real power aggregate of Node ", num2str(i)))
%     
%     ax3 = nexttile;
%     plot(ax3,q_aggregate{i})
%     hold on
%     plot(ax3,zeros(48,1),'k--')
%     grid on
%     title(ax3, strcat("Reactive power aggregate of Node ", num2str(i)))
% 
% end

%% CADMM parameters
% Censoring parameters
censoring_rho = 0.98;
censoring_alpha = 1e-12;
censoring_r = 1e-12;

%% Censored ADMM Solution

% initialise the variables
xi_cens = zeros(2*K*T, N);
u_cens = zeros(2*T+2*K*T,N);
lambda_cens = zeros(2*K*T,N);
nu_cens = zeros(2*K*T,N);

u_state_cens = zeros(2*T+2*K*T,N);
lambda_state_cens = zeros(2*K*T,N);
nu_state_cens = zeros(2*K*T,N);

% initialise storage
lambda_cens_step = {lambda_cens};
nu_cens_step = {nu_cens};
u_cens_step = {u_cens};
lambda_state_cens_step = {lambda_state_cens};
nu_state_cens_step = {nu_state_cens};
u_state_cens_step = {u_state_cens};


cens_iter = 0;
change_cens = 1;
U_cens_err = [];
Lambda_cens_err = [];
Nu_cens_err = [];
Change_cens_step = [change_cens];
time_cens = [];

xi_step = {};

Transmission = zeros(N, iteration_limit);

disp("Computing censored solutions")

while ((change_cens>progress_tol) || isnan(change_cens)) && (cens_iter<iteration_limit) 
    % Update u and lambda
    cens_iter = cens_iter+1;
    tic;

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
            temp = temp + ( lambda_state_cens(:,i) + lambda_state_cens(:,neighbors{i}(j)) );
        end
        
        obj = (eta'*p + kappa*(p'*p)) +...
              (step_size/(n_deg(i)*4))*norm( (1/step_size)*(Psi{i}*(x) - (1/N)*w) - (1/step_size)*nu_cens_old(:,i) +  temp)^2;
        
        options = sdpsettings('verbose',0);
        optimize(constraints, obj, options);
        u_cens(:,i) = value(x); 
        
        lambda_cens(:,i) = (1/(2*n_deg(i)))*( temp - (1/step_size)*nu_cens_old(:,i) + (1/step_size)*( Psi{i}*u_cens(:,i) - (1/N)*w )  ); 
    end 
    
    u_cens_step{cens_iter} = u_cens;
    lambda_cens_step{cens_iter} = lambda_cens;
    
    % Transmission Step
    iter_transmissions = 0;
    for i=1:N
        % Computing the difference between current state and primal update
        xi(:,i) = lambda_state_cens(:,i) - lambda_cens(:,i);
%         disp(norm(xi(:,i))^2)
        % Transmission loop
        if ( norm(xi(:,i))^2 - censoring_alpha*( censoring_rho^cens_iter ) )>=0
            lambda_state_cens(:,i) = lambda_cens(:,i);
            Transmission(i,cens_iter) = 1;
            iter_transmissions = iter_transmissions + 1;
        end
    end
    
%     xi_step{cens_iter} = xi;
    
    lambda_state_cens_step{cens_iter} = lambda_state_cens;
    
    % Updating nu
    parfor i=1:N
        temp1 = 0;
        for j=1:length(neighbors{i})
            temp1 = temp1 + lambda_state_cens(:,i) - lambda_state_cens(:,j);
        end
        
        nu_cens(:,i) = nu_cens_old(:,i) + step_size*temp1;
    end

    nu_cens_step{cens_iter} = nu_cens;
    
    u_cens_err = norm(u_cens_old-u_cens,'fro')^2;
    nu_cens_err = norm(nu_cens_old-nu_cens,'fro')^2;
    lambda_cens_err = norm(lambda_cens_old - lambda_cens, 'fro')^2;
    
    U_cens_err = [U_cens_err, u_cens_err];
    Lambda_cens_err = [Lambda_cens_err, lambda_cens_err];
    Nu_cens_err = [Nu_cens_err, nu_cens_err];
    
    change_cens=max([u_cens_err, nu_cens_err, lambda_cens_err]);
    Change_cens_step = [Change_cens_step, change_cens];
    
    time_taken = toc;
    time_cens = [time_cens, time_taken];
    
    disp(strcat("Citer - ",int2str(cens_iter)," | prog - ", num2str(round(change_cens,11)), " | # ", num2str(iter_transmissions), " | time - ", num2str(round(time_taken/60,1))))
    
end

beep;

%% Plotting Solutions

% disp(V_baseline)

V_soln = zeros(K*T,1);

p_ADMM = {};
q_ADMM = {};
p_cens = {};
q_cens = {};
for i=1:N
    p_ADMM{i} = u_ADMM(   1:T   , i);
    q_ADMM{i} = u_ADMM( T+1:2*T , i);
    p_cens{i} = u_cens(   1:T   , i);
    q_cens{i} = u_cens( T+1:2*T , i);
end

V_control_ADMM = zeros(K*T,1);
V_control_cens = zeros(K*T,1);

for i=1:N
    V_control_ADMM = V_control_ADMM + D_stack{i}*p_ADMM{i} + E_stack{i}*q_ADMM{i};
    V_control_cens = V_control_cens + D_stack{i}*p_cens{i} + E_stack{i}*q_cens{i};
end

V_soln_ADMM = V_stacked + (V_control_ADMM/(v0^2));
V_soln_cens = V_stacked + (V_control_cens/(v0^2));

V_ADMM_reshaped = reshape(V_soln_ADMM, [K,T]);
V_cens_reshaped = reshape(V_soln_cens, [K,T]);

p_aggregate = {};
q_aggregate = {};

for i=1:K
    p = zeros(48,1);
    q = zeros(48,1);
    if i==1
        for j=1:30
            p = p + p_ADMM{j};
            q = q + q_ADMM{j};
        end
    elseif i==2
        for j=31:60
            p = p + p_ADMM{j};
            q = q + q_ADMM{j};
        end
    elseif i==3
        for j=61:90
            p = p + p_ADMM{j};
            q = q + q_ADMM{j};
        end
    end
    p_aggregate{i} = p;
    q_aggregate{i} = q;
end

figure()
plot(sqrt(V_b(1,:)), '-','color',[0.4660 0.6740 0.1880], 'DisplayName', 'Baseline Voltage - Phase 1')
hold on
plot(sqrt(V_ADMM_reshaped(1,:)),'-^','color',[0.4660 0.6740 0.1880], 'DisplayName', 'Dis-Net-EVCD - Phase 1')
hold on
plot(sqrt(V_cens_reshaped(1,:)),'-p','color',[0.4660 0.6740 0.1880], 'DisplayName', 'Algorithm 1 - Phase 1')
hold on
plot(sqrt(V_b(2,:)), '-r', 'DisplayName', 'Baseline Voltage - Phase 2')
hold on
plot(sqrt(V_ADMM_reshaped(2,:)),'-r^', 'DisplayName', 'Dis-Net-EVCD - Phase 2')
hold on
plot(sqrt(V_cens_reshaped(2,:)),'-rp', 'DisplayName', 'Algorithm 1 - Phase 2')
hold on
plot(sqrt(V_b(3,:)), 'b-', 'DisplayName', 'Baseline Voltage - Phase 3')
hold on
plot(sqrt(V_ADMM_reshaped(3,:)),'b^-', 'DisplayName', 'Dis-Net-EVCD - Phase 3')
hold on
plot(sqrt(V_cens_reshaped(3,:)),'bp-', 'DisplayName', 'Algorithm 1 - Phase 3')
grid on
plot((1+percentage)*ones(1,T), 'k--')
hold on
plot((1-percentage)*ones(1,T), 'k--')
hold on
title("Voltage Profile at Node 1")
lgd = legend('Baseline Voltage - Phase 1', 'Benchmark - Phase 1', 'CC-ADMM - Phase 1', ...
    'Baseline Voltage - Phase 2', 'Benchmark - Phase 2', 'CC-ADMM - Phase 2', ...
    'Baseline Voltage - Phase 3', 'Benchmark - Phase 3', 'CC-ADMM - Phase 3', '', '');
lgd.NumColumns=3;
lgd.FontSize = 20;
ylabel("Nodal Voltage (p.u.)")
xlim([1,48])
xlabel("Time of day")
xticks([1, 12, 24, 36, 48])
xticklabels(["12am", "6am", "12pm", "6pm", "12am"])
yticks([0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06])
yticklabels([" ", "0.95"," "," "," "," ","1.00"," "," "," "," ","1.05"," "])
set(findall(gcf,'-property','FontSize'),'FontSize',20)

%     tiledlayout(3,1)
    
%     ax1 = nexttile;
%     plot(ax1,sqrt(V_b(i,:)), 'rx--')
%     hold on
%     plot(ax1,(1+percentage)*ones(1,T), 'k-')
%     hold on
%     plot(ax1,(1-percentage)*ones(1,T), 'k-')
%     hold on
%     plot(ax1,sqrt(V_ADMM_reshaped(i,:)),'g--')
%     hold on
%     plot(ax1,sqrt(V_cens_reshaped(i,:)),'bo')
%     grid on
%     title(ax1, strcat("Baseline Voltage and Control Voltage of Node ", num2str(i)))
%     
%     ax2 = nexttile;
%     plot(ax2,p_aggregate{i})
%     hold on
%     plot(ax2,zeros(48,1),'k--')
%     grid on
%     title(ax2, strcat("Real power aggregate of Node ", num2str(i)))
%     
%     ax3 = nexttile;
%     plot(ax3,q_aggregate{i})
%     hold on
%     plot(ax3,zeros(48,1),'k--')
%     grid on
%     title(ax3, strcat("Reactive power aggregate of Node ", num2str(i)))

% end

%% Error Plots
% 
% col = parula(iter); % rgb values of 10 colors
% 
% for i=1:3
%     figure()
%     for k=1:iter
%         plot(lambda_ADMM_step{k}(:,i), 'color',col(k,:))
%         hold on
%         grid on
%     end
%     legend('location','northwest')
% end
% 
% col = parula(cens_iter);
% for i=1:3
%     figure()
%     for k=1:cens_iter
%         plot(lambda_cens_step{k}(:,i), 'color',col(k,:))
%         hold on
%         grid on
%     end
%     legend('location','northwest')
% end
% 
% % Error_U = [];
% % Error_L = {};
% % Error_N = [];
% % Error_LS = [];
% % 
% % 
% % for j=1:N
% %     Error_L{j} = [];
% %     for i=1:cens_iter-1
% %         
% %         %     Error_U = [Error_U, norm( u_cens_step{i} - u_cens_step{i+1}, 'fro' )^2];
% %         Error_L{j} = [Error_L{j}, norm( lambda_cens_step{i}(:,j) - lambda_cens_step{i+1}(:,j))^2];
% %         %     Error_N = [Error_N, norm( nu_cens_step{i} - nu_cens_step{i+1}, 'fro' )^2];
% %         %     Error_LS = [Error_LS, norm( lambda_state_cens_step{i}-lambda_state_cens_step{i+1} )^2];
% %     end
% % end
% % 
% % cens_curve = [];
% % 
% % for i=1:cens_iter
% %     cens_curve = [cens_curve, censoring_alpha*(censoring_rho^i)];
% % end
% % 
% % 
% % figure()
% % % plot(Error_U, 'rx-')
% % % hold on
% % for j=1:N
% %     plot(Error_L{j}, 'bx-')
% %     hold on 
% % % plot(Error_N, 'kx-')
% % % hold on
% % % plot(Error_LS, 'bo')
% % % % hold on
% % %     plot(cens_curve, 'g-')
% % %     legend('lambda', 'cencoring')
% % end
% % 
% % % lambda_1 = [];
% % % for i=1:cens_iter
% % %     lambda_1 = [lambda_1, lambda_cens_step{i}(144,1)];
% % % end
% % % figure()
% % % plot(lambda_1)

%% sum of transmissions

sum_trans = 0;
for i=1:N
    sum_trans = sum_trans + sum(Transmission(i,:));
end

figure()
clims=[0,1];
imagesc(Transmission, clims);
colormap(gray(256));
title("Communication Censoring Pattern")
xlabel("Iteration")
ylabel("Customer ID")
yticks([0:50:150])
xticks([0:10:30])
set(findall(gcf,'-property','FontSize'),'FontSize',30)