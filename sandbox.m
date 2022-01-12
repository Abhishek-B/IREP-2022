%{
Censored ADMM code in the primal setting from the paper
"COMMUNICATION-CENSORED ADMM FOR DECENTRALIZED CONSENSUS OPTIMIZATION"
-Abhishek B. Dec. 2021
%}

set(0,'DefaultFigureWindowStyle','docked')

%% Graph parameters setup

n = 10; % number of nodes/customers whatever
p = 3; % dimension of primal variable for each node

Adjacency = ones(n,n) - eye(n,n);


%% cost function
% working for example A of the paper, which is least squares cost function

% choosing the true solution
xtrue = rand(p,1);
Xtrue = [];
for i=1:n
    Xtrue = horzcat(Xtrue, xtrue);
end

% choosing the random data for each node.
% for test, will set A to identity and y to be the true solution

A = {};
y = {};
f = {};

for i=1:n
    A{i} = eye(p,p); %10*rand(p,p); %
    y{i} = xtrue;
    f{i} = @(x) (A{i}*x - y{i})'*(A{i}*x-y{i});
end


%% Classic ADMM Solve
iterations = 250;
step_size = 2^(0);
[xadmm,X,zadmm,Z,ladmm,L] = ADMM(iterations,n,p,step_size,f);

%% Censored ADMM Solve
% Will use same iteration count as ADMM, but will vary the stepsize to see
% changes

step_size_c = 0.35;
rho=0.5;
alpha=0.1; % with alpha=0 this should reduce to simple decentralized admm which should converge...

[xc,Xc, x_state,xi,transmission,lambda,Lc] = ADMM_censored(iterations,n,p,rho,alpha,step_size_c,Adjacency,f);


%% Error comparison

err=[];
errc = [];
err_iter=[];
for i=1:iterations
    err = [err, ( norm(X{i}-Xtrue,'fro')^2 / norm(Xtrue,'fro')^2 )];
    errc = [errc, ( norm(Xc{i}-Xtrue,'fro')^2 / norm(Xtrue,'fro')^2 )];
%     err_iter = [err_iter, ( norm(x_iter{i}-Xtrue,'fro')^2 / norm(Xtrue,'fro')^2 )];
end

semilogy([1:iterations],err,'k')
grid on
hold on
semilogy([1:iterations],errc,'r')
% hold on
% semilogy([1:iterations],err_iter,'bx')

%% Implementing classic consensus ADMM algorithm of Boyd to compare things to

% x_ADMM = zeros(p,n);
% l_ADMM = zeros(p,n);
% z_ADMM = zeros(p,1);
% rho_ADMM = 1;
% 
% C = {};
% for i=1:n
%     C{i} = A{i}'*A{i}+rho_ADMM*eye(p,p);
% end
% 
% for k=1:iterations
%     for i=1:n
%         x_ADMM(:,i) = C{i}\( A{i}*y{i} + z_ADMM - 0.5*l_ADMM(:,i) );
%     end
% 
% 
%     temp3 = 0;
%     for i=1:n
%         temp3 = temp3+ x_ADMM(:,i) + (1/rho_ADMM)*l_ADMM(:,i);
%     end
%     z_ADMM = (1/n)*temp3;
% 
%     for i=1:n
%         l_ADMM(:,i) = l_ADMM(:,i) + rho_ADMM*(x_ADMM(:,i) - z_ADMM);
%     end
% end


%% Update step for censored ADMM
% for this cost function the primal update has a closed form expression, so
% we will use that here. In other examples this would be replaced with a
% fminsearch thingy

% G = graph(Adjacency);
% 
% edge = G.Edges;
% 
% n_deg = degree(G);
% 
% neighbors = {};
% for i=1:n
%     neighbors{i} = find(Adjacency(:,i));
% end
% 
% c = 20; % stepsize parameter of COCA (should be tuned later)
% alpha = 0.1; % parameter for determining transmission
% rho = 2^(-10); % sequence for determining transmission
% 
% x = zeros(p,n); % primal value for each customer
% x_state = zeros(p,n); % state value for each customer, the neighbors states
%                       % can be found from different columns of this.
%                       
% x_iter={};
% 
% xi = zeros(p,n); % the difference between primal value and state value
% lambda = zeros(p,n); % dual variabels for each customer
% 
% B = {};
% for i=1:n
%     B{i} = (A{i}'*A{i} + 2*c*n_deg(i)*eye(p,p)); % Not using inv() since that is not as good
% end
% 
% transmission = zeros(n,iterations);
% 
% err2 = [];
% 
% for k=1:iterations
%     
%     for i=1:n
%         temp=0;
%         for j=1:length(neighbors{i})
%             temp = temp +( x_state(:,i) - x_state(:,neighbors{i}(j)) );
%         end
%         
%         x(:,i) = B{i}\( A{i}'*y{i} - lambda(i) + c*temp );
%         xi(:,i) = x_state(:,i) - x(:,i); 
%         
%         if ( norm(xi(:,i))^2 - alpha*(rho^k) )>=0
%             %disp(norm(xi(:,i))^2 - alpha*(rho^k))
%             %disp("updating transmission")
%             x_state(:,i) = x(:,i);
%             transmission(i,k) = 1;
%         end
%     end
%     
%     err2 = [err2, ( norm(x - Xtrue,'fro')^2 / norm(Xtrue,'fro')^2   )];
%     
%     for i=1:n
%         temp2 = 0;
%         for j=1:length(neighbors{i})
%             temp2 = temp2 + ( x_state(:,i) - x_state(:,neighbors{i}(j)) );
%         end
%         lambda(:,i) = lambda(:,i) + c*temp2;
%         %disp(lambda(:,i))
%     end
%     
%     x_iter{k} = x;
% end

% Stuff is working now, but with all the parameter tuning...
% who knows what's going on



