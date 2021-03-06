function [x,X,x_state,xi,transmission,lambda,L] = ADMM_censored(iteration,n,p,rho,alpha,c,Adjacency,costs)
%UNTITLED Summary of this function goes here
% c = 20; % stepsize parameter of COCA (should be tuned later)
% alpha = 0.1; % parameter for determining transmission
% rho = 2^(-10); % sequence for determining transmission



%% Graph Parameters

G = graph(Adjacency);

edge = G.Edges;

n_deg = degree(G);

neighbors = {};
for i=1:n
    neighbors{i} = find(Adjacency(:,i));
end
%%
x = zeros(p,n);       % primal value for each customer
x_state = zeros(p,n); % state value for each customer, the neighbors states
                      % can be found from different columns of this.

X = {}; % cell to hold each iteration value of the primal
L = {}; % cell to hold each iteration value of the dual

xi = zeros(p,n);      % the difference between primal value and state value
lambda = zeros(p,n);  % dual variabels for each customer

transmission = zeros(n,iteration);

for k=1:iteration
    
    % Update of the primal variable
    for i=1:n
        % Building the sum of difference between self state and neighbors
        % states
        temp=0;
        for j=1:length(neighbors{i})
            temp = temp +( x_state(:,i) + x_state(:,neighbors{i}(j)) );
        end
        
        % the overall cost function for the primal update
        fun = @(x)costs{i}(x) + x'*(lambda(:,i) - c*temp) + c*n_deg(i)*norm(x)^2;
        
        % updating the primal through a subminimization
        x(:,i) = fminsearch(fun, zeros(p,1));
        
        % Computing the difference between current state and primal update
        xi(:,i) = x_state(:,i) - x(:,i); 
        %disp(xi(:,i)')
        % Transmission loop
        if ( norm(xi(:,i))^2 - alpha*(rho^k) )>=0
            %disp(norm(xi(:,i))^2 - alpha*(rho^k))
            %disp("true")
            x_state(:,i) = x(:,i);
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
            temp2 = temp2 + ( x_state(:,i) - x_state(:,neighbors{i}(j)) );
        end
        
        % Updating the dual variable
        lambda(:,i) = lambda(:,i) + c*temp2;
    end
    
    % Storing the primal and dual updates for the iteration
    X{k} = x;
    L{k} = lambda;
end





























end

