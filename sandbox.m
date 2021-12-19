%{
Censored ADMM code in the primal setting from the paper
"COMMUNICATION-CENSORED ADMM FOR DECENTRALIZED CONSENSUS OPTIMIZATION"
-Abhishek B. Dec. 2021
%}

%% Graph parameters setup
% will need to fix all this later into a functions of sorts, for now will
% just write down for a simple example.

n = 3; % number of nodes/customers whatever
p = 3; % dimension of primal variable for each node
c = 1; % stepsize parameter of COCA (should be tuned later)
alpha = 1; % parameter for determining transmission
rho = 2^(-10); % sequence for determining transmission

Adjacency = ones(n,n) - eye(n,n);
G = graph(Adjacency);

edge = G.Edges;

n_deg = degree(G);

neighbors = {};
for i=1:n
    neighbors{i} = find(Adjacency(:,i));
end

%%

x = zeros(n,p); % primal value for each customer
x_state = zeros(n,p); % state value for each customer, the neighbors states
                      % can be found from different columns of this.

xi = zeros(n,p); % the difference between primal value and state value
lambda = zeros(n,p); % dual variabels for each customer


%% cost function
% working for example A of the paper, which is least squares cost function

A = {};
y = {};
f = {};

for i=1:3
    A{i} = eye(p,p); %10*rand(p,p);
    y{i} = ones(p,1); %10*rand(p,1);
    f{i} = @(x) (A{i}*x - y{i})'*(A{i}*x-y{i});
end
%% Implementing classic consensus ADMM algorithm of Boyd to compare things to

iterations = 1000;

x_ADMM = zeros(n,p);
l_ADMM = zeros(n,p);
z_ADMM = zeros(p,1);
rho_ADMM = 2^(-10);

C = {};
for i=1:n
    C{i} = A{i}'*A{i}+rho_ADMM*eye(p,p);
end

for k=1:iterations
    for i=1:n
        x_ADMM(:,i) = C{i}\( A{i}*y{i} + z_ADMM - 0.5*l_ADMM(:,i) );
    end
end

temp3 = 0;
for i=1:n
    temp3 = temp3+ x_ADMM(:,i) + (1/rho_ADMM)*l_ADMM(:,i);
end
z_ADMM = (1/n)*temp3;

for k=1:iterations
    for i=1:n
        l_ADMM(:,i) = l_ADMM(:,i) + rho_ADMM*(x_ADMM(:,i) - z_ADMM);
    end
end

%% Update step for censored ADMM
% for this cost function the primal update has a closed form expression, so
% we will use that here. In other examples this would be replaced with a
% fminsearch thingy

B = {};
for i=1:n
    B{i} = (A{i}'*A{i} + 2*c*n_deg(i)*eye(p,p)); % Not using inv() since that is not as good
end

%%


transmission = zeros(n,iterations);

for k=1:iterations
    for i=1:n
        temp=x_state(:,i);
        for j=1:length(neighbors{i})
            temp = temp - x_state(:,neighbors{i}(j));
        end
        
        x(:,i) = B{i}\( A{i}'*y{i} - lambda(i) + c*temp );
        
        %disp(x(:,i))
        
        xi(:,i) = x_state(:,i) - x(:,i); 
        
        %disp(xi(:,i))
        
        %disp(norm(xi(:,i))^2 - alpha*(rho^k))
        
        if ( norm(xi(:,i))^2 - alpha*(rho^k) )>=0
            %disp(norm(xi(:,i))^2 - alpha*(rho^k))
            %disp("updating transmission")
            x_state(:,i) = x(:,i);
            transmission(i,k) = 1;
        end
    end
    
    for i=1:n
        temp2 = x_state(:,i);
        for j=1:length(neighbors{i})
            temp2 = temp2 - x_state(:,neighbors{i}(j));
        end
        lambda(:,i) = lambda(:,i) + c*temp2;
        %disp(lambda(:,i))
    end
end




