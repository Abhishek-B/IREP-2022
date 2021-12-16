%{
Censored ADMM code in the primal setting from the paper
"COMMUNICATION-CENSORED ADMM FOR DECENTRALIZED CONSENSUS OPTIMIZATION"
-Abhishek B. Dec. 2021
%}

%% Graph parameters setup
% will need to fix all this later into a functions of sorts, for now will
% just write down for a simple example.

nodes = 3;
edges = {[1,2], [1,3], [2,3]};

%% COCA algorithm run by node i

x = zeros(nodes,1);
x_state = zeros(nodes,1);
l = zeros(nodes,1);
for i=1:nodes
    for k=1:iterations
        x(i) = primal_update();
        xi = x_state(i) - x(i);
        if H(k,xi)>=0
            x_state(i) = x(i);
            %Trannsmit to neighbors
        end
        