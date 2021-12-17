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
rho = 0.5; % sequence for determining transmission

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
    A{i} = 10*rand(p,p);
    y{i} = 10*rand(p,1);
    f{i} = @(x) (A{i}*x - y{i})'*(A{i}*x-y{i});
end


%% Update step
% for this cost function the primal update has a closed form expression, so
% we will use that here. In other examples this would be replaced with a
% fminsearch thingy

B = {};
for i=1:n
    B{i} = inv((A{i}'*A{i} + 2*c*n_deg(i)*eye(p,p)));
end

iterations = 1000;

for k=1:iterations
    transmission = zeros(n,1);
    for i=1:n
        temp=x_state(:,i);
        for j=1:length(neighbors{i})
            temp = temp - x_state(:,neighbors{i}(j));
        end
        
        x(:,i) = B{i}*( A{i}'*y{i} - lambda(i) + c*temp );
        
        xi(:,i) = x_state(:,i) - x(:,i); 
        
        if (norm(xi(:,i)) - alpha*(rho^k))>=0
            x_state(:,i) = x(:,i);
        end
        transmission(i) = 1;
    end
    
    for i=1:n
        temp2 = x_state(:,i);
        for j=1:length(neighbors{i})
            temp2 = temp2 - x_state(:,neighbors{i}(j));
        end
        lambda(:,i) = lambda(:,i) + c*temp2;
    end
end