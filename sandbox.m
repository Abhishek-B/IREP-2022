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


%% Primal update step
% for this cost function the primal update has a closed form expression, so
% we will use that here. In other examples this would be replaced with a
% fminsearch thingy











%% COCA algorithm run by node i

x = zeros(n,1);
x_state = zeros(n,1);
l = zeros(n,1);
for i=1:n
    for k=1:iterations
        x(i) = primal_update();
        xi = x_state(i) - x(i);
        if H(k,xi)>=0
            x_state(i) = x(i);
        end
    end
end


% Transmission Step


% Dual variable update