%{
Testing the IEEE 13 node feeder.
Jan 2022 - Abhishek B.
%}
set(0,'DefaultFigureWindowStyle','docked')
%% REES model parameters and such

% Customers, time steps, supply point numbers
N = 23; % 600 is too much for the purpose of demonstration, there are 29 
        % supply points with 6 having 0 customers, so I will put 1 customer
        % on each supply point
T = 48;
K = 29;

% power profile vectors (empty for now)
% maybe I dont need these here...eh whatever, I'll come back to it
% p = zeros(T,N);
% q = zeros(T,N);
% s = zeros(2*K*T, N);
% u = [p;q;s];

% Stacking design parameters into row vectors, each elements corresponds to
% customer i
% setting to zero for now, will fill in soon

% Changing these to simplest possible scenario.
a = 0; %randi([10,14],1,N); % arrival time 5pm-7pm
d = T; %randi([38,42],1,N); % departure time 7am-9am

capacities = [42,21,24,20,30,16,23,16.5,28,60,90,75];
b = zeros(1,N); 
for i=1:N
    b(i) = capacities(randi([1,12]));
end
s_bar = 10; % same for all customers %*ones(1,N); % inverter capacity
sig_hat = (0.3+0.2*rand(1,N)).*b; % initial SoC
sig_u = 0.85*b; % max SoC
sig_l = 0.2*b;  % min SoC
sig_star = zeros(1,N); % target SoC - this isnt listed in the paper, so I 
                       % will set this to between sig_hat and sig_upper,
                       % rounded to the nearest integer
for i=1:N
    sig_star(i) = ceil(sig_hat(i) + (sig_u(i)-sig_hat(i))*rand);
end
p_u = 6.6;  % max charge rate
p_l = -6.6; % min charge rate

% One of the places where these get used, Nanduni has missed, so I will set
% both to 1, to see if Yalmip will work.
mu_u = 1; %0.9; % charge efficiency
mu_l = 1; %1.1; % discharge efficiency


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

% Below comment was for 600 customers, I am changing to now 23.
% % The only way I can make sense of the numbers in Nanduni's paper is if
% % each supply point has its own number of customers. So each customer is
% % essentially only connected to a single phase...
% % yeah this incrementation and labelling of customers is a bit confusing...

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

% theta(4,1:26) = 1;
% theta(5,27:52) = 1;
% theta(6,53:78) = 1;
% theta(7,79:104) = 1;
% theta(8,105:130) = 1;
% theta(9,131:156) = 1;
% theta(10,157:182) = 1;
% theta(11, 183:208) = 1;
% theta(12, 209:234) = 1;
% theta(13, 235:260) = 1;
% theta(17, 261:286) = 1;
% theta(18, 287:312) = 1;
% theta(19, 313:338) = 1;
% theta(20, 339:364) = 1;
% theta(21, 365:390) = 1;
% theta(22, 391:416) = 1;
% theta(23, 417:442) = 1;
% theta(24, 443:468) = 1;
% theta(25, 469:495) = 1; % Supply point 10 which has 27 customers
% theta(26, 496:522) = 1; % Supply point 11 which has 27 customers
% theta(27, 523:548) = 1;
% theta(28, 549:574) = 1;
% theta(29, 575:600) = 1;

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


%% Old remnant code for baseline voltages
 
% baseline voltage
% I'm gonna set this up myself by choosing random values. The dataset from
% ausgrid doesn't seem to have the link up anymore, and i wouldn't know how
% to convert nodal voltage data to the p.u. system. 


% these have 49 values and not 48...so maybe the time steps are 49=T
x = [   1,    4,     5, 8,   11,    16,   17,    20,   23,   32,    34,    37,   40, 42,   45, 48];
xq = [1:1:48];
v = [1.04, 1.01, 1.023, 1, 0.95, 0.955, 0.96, 0.974, 0.98, 0.98, 0.976, 0.974, 0.97,  1, 1.01,  1];

v_temp = interp1(x,v,xq);

% figure()
% plot(xq, v_temp, 'r-o')
% ylim([0.9,1.1])

V_baseline = zeros(K,T);

for j=1:T
    for i=1:K
        V_baseline(i,j) = [v_temp(j)-0.01*v_temp(j)] + 0.02*rand;
    end
end

% figure()
% for i=1:29
%     plot(xq, V_baseline(i,:), 'k-o')
%     ylim([0.9,1.1])
%     hold on
%     grid on
%     set(gca,'ytick',[0.0:0.01:1.1])
% end
% plot(xq, 1.046*ones(1,48), 'r-')
% plot(xq, 0.954*ones(1,48), 'r-')

% The plot seem okay-ish for these baseline voltages. Eyeballing it, they
% seem roughly similar to Nanduni's plot, so I'll go with this.

% V_base = [];
% for i=1:T
%     V_base = [V_base, V_baseline(:,i)];
% end

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
%      -v_lower + V_base_stacked];

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

% costs={};
% for i=1:N
%     costs{i} = @(u) eta'*u(1:T) + kappa*u(1:T)'*u(1:T); %Omega_n(u, T, eta, kappa); % + Indicator_n(u, T, w, s_bar(i), mu_u, mu_l, a(i), d(i), p_u, p_l, alpha_u(i), alpha_l(i), e(i));
% end

% Do i really need this crap if I'm using YALMIP?...
% nah I dont think I do

%% ADMM
tic;
iteration = 10;
lambda = zeros(2*K*T,N);       % primal value for each customer
lambda_state = zeros(2*K*T,N); % state value for each customer, the neighbors states
                      % can be found from different columns of this.

c_stepsize = 0.0002;
alpha = 0;
rho=0;
lambda_step = {}; % cell to hold each iteration value of the primal
nu_step = {}; % cell to hold each iteration value of the dual

xi = zeros(2*K*T,N);      % the difference between primal value and state value
nu = zeros(2*K*T,N);  % dual variabels for each customer
U = zeros(2*(T+K*T),N);

transmission = zeros(N,iteration);


I = eye(T,T);
O = ones(T,1);
Z = zeros(T,1);

% Building vector mu, of charge efficiencies
mu = ones(T,1); %zeros(T,1)
%         for j=1:T
%             if p(j)>=0
%                 mu(j)=mu_u;
%             else
%                 mu(j)=mu_l;
%             end
%         end
% Building matrix M which bounds charge rates
M = tril(ones(T,T));
% for j=1:T
%     for l=1:T
%         if j>=l
%             M(j,l) = mu(l);
%         else
%             M(j,l)=0;
%         end
%     end
% end

for k=1:iteration
    disp(k)
    % looping over customers
    parfor i=1:N
        
%         disp([k, i])
        
        % Building the sum of self state and neighbors
        % states
        temp=0;
        for j=1:length(neighbors{i})
            temp = temp +( lambda_state(:,i) + lambda_state(:,neighbors{i}(j)) );
        end
        
        %% Updating u_n for the primal step
        % Using YALMIP for this one
        
        %% Setup of SDPVARS
        % each customers u_n
        x = sdpvar(2*(T+K*T),1);
        
        % finding p,q,s from u
        p = x(1:T);
        q = x(T+1:2*T);
        s = x(2*T+1:end);
%         
%         % sanity check
%         if ~size(p,1)==size(q,1)
%             error("sizes don't match for p and q")
%         end
        
        %% constraints - inverter capacity
        
        constraints = [p.^2 + q.^2 <= s_bar*ones(T,1)];
        %inv_cap = s_bar(i)^2;
%         for j=1:T
%             constraints = [constraints, p(j)^2+q(j)^2 <= s_bar]; %inv_cap];
%         end
        
        %% Constraints building other stuff

        constraints = [constraints, M*p <= alpha_u(i)*O, alpha_l(i)*O <= M*p];
        
        constraints = [constraints, p_l*O <= p, p<= p_u*O];
        
        constraints = [constraints, O'*p >= e(i)];
        
        % this part is not needed when I set a=0, d=T.
%         % Building availibility matrix L
%         L = zeros(T,T);
%         for j=1:T
%             for l=1:T
%                 if (j==l)&&(a(i)<j)&&(j<=d(i))
%                     L(j,l)=1;
%                 end
%             end
%         end
        

%         F = [I ; -I ; M ; -M ; mu' ; -mu' ; (I-L) ; (L-I)];
%         f = [p_l*O ; -p_u*O ; alpha_l(i)*O ; -alpha_u(i)*O ; e(i) ; -e(i) ; Z ; Z];
         
        %% Constraints - final constraint for feasibility set
        % Final feasibility set constraint for (p,q), and then also for the
        % slack variable 0<=s<=w
%         constraints = [constraints, f<=F*p, zeros(length(s),1)<=s, s<=w];
        
        constraints = [constraints, zeros(length(s),1) <= s]; %, s<=w];
        %% Updating u_n with optimize and Mosek
        
        obj = eta'*p + kappa*(p'*p);
        extr_obj = (1/(4*c_stepsize*n_deg(i)))*norm( Psi{i}*x - (1/N)*w - nu(:,i) + c_stepsize*temp )^2  ; 
        opt = sdpsettings('verbose',0);
        optimize(constraints, obj+extr_obj, opt);
        u = value(x);
        U(:,i) = u;
        
        
        %% Updating lambda
        
        
        lambda(:,i) = (1/(2*c_stepsize*n_deg(i)))*( Psi{i}*u - (1/N)*w - nu(:,i) + c_stepsize*temp );
        
%         
%         % the overall cost function for the primal update
%         fun = @(x)costs{i}(x) + x'*(lambda(:,i) - c*temp) + c*n_deg(i)*norm(x)^2;
%         %%
%         % updating the primal through a subminimization
%         lambda(:,i) = fminsearch(fun, zeros(p,1));
%         
        
%         %% Transmission Step
%         % Computing the difference between current state and primal update
%         xi(:,i) = lambda_state(:,i) - lambda(:,i); 
%         % Transmission loop
%         if ( norm(xi(:,i))^2 - alpha*(rho^k) )>=0
%             lambda_state(:,i) = lambda(:,i);
%             transmission(i,k) = 1;
%         end
    end
    
    %% Dual variable update
    parfor i=1:N
        
        % building the sum of the differences between self state and
        % neighbors current states
        temp2 = 0;
        for j=1:length(neighbors{i})
            temp2 = temp2 + ( lambda_state(:,i) - lambda_state(:,neighbors{i}(j)) );
        end
        
        % Updating the dual variable
        nu(:,i) = nu(:,i) + c_stepsize*temp2;
    end
    
    %% Storing the primal and dual updates for the iteration
    lambda_step{k} = lambda;
    nu_step{k} = nu;
end

time_taken=toc;
%% Reconstruct Voltages from solutions

V_soln = zeros(K*T,1);

temp = zeros(K*T,1);
for i=1:N
    temp = temp + D_stack{i}*U(1:T, i) + E_stack{i}*U(T+1:2*T, i);
end

V_soln = sqrt(V_baseline(:) + temp);

V_soln_reshaped = reshape(V_soln, [K,T]);

% V_soln = zeros(K,T);
% padding = zeros(T,1);
% p_soln = [padding, padding, padding, U(1:T,1:10), padding, padding, padding, U(1:T,11:end)];
% q_soln = [padding, padding, padding, U(T+1:2*T,1:10), padding, padding, padding, U(T+1:2*T,11:end)];
% 
% 
% for i=1:T
%     V_soln(:,i) = V_base(:,i) - R*p_soln(i,:)' - X*q_soln(i,:)';
% end

%% Plots

% figure()
% for i=1:29
%     plot(xq, V_baseline(i,:), 'k-o')
%     ylim([0.9,1.1])
%     hold on
%     grid on
%     set(gca,'ytick',[0.0:0.01:1.1])
% end
% plot(xq, 1.046*ones(1,48), 'r-')
% plot(xq, 0.954*ones(1,48), 'r-')
figure()
for i=1:K
    plot(V_baseline(i,:), 'ro--')
    ylim([0.9,1.1])
    hold on
    grid on
    set(gca,'ytick',[0.0:0.01:1.1])
    hold on
    plot(V_soln_reshaped(i,:), 'go-')
end
plot(1.046*ones(1,48), 'k-')
plot(0.954*ones(1,48), 'k-')
