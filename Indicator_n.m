function [y] = Indicator_n(p,q,s,w,Psi, s_bar,mu_u,mu_l, a,d, p_u,p_l, alpha_u,alpha_l,e)
%Indicator function for the battery constraints

% the constants alpha_u,alpha_l,e can be computed 
% in script

if ~size(p,1)==size(q,1)
    error("sizes don't match for p and q")
end

N = size(p,2);
T = size(p,1);

I = eye(T,T);

% Checking apparent power constraint.
for i=1:T
    if ~(p(i)^2+q(i)^2<=s_bar^2)
        y=Inf;
        return;
    end
end

% Building vector mu, of charge efficiencies
mu = zeros(T,1);
for i=1:T
    if p(i)>=0
        mu(i)=mu_u;
    else
        mu(i)=mu_l;
    end
end

% Building matrix M which bounds charge rates
M = zeros(T,T);
for i=1:T
    for j=1:T
        if i>=j
            M(i,j) = mu(j);
        else
            M(i,j)=0;
        end
    end
end

% Building availibility matrix L
L = zeros(T,T);
for i=1:T
    for j=1:T
        if (i==j)&&(a<i<=d)
            L(i,j)=1;
        end
    end
end

%% building final matrices

O = ones(T,1);
Z = zeros(T,1);

F = [I ; -I ; M ; -M ; mu' ; -mu' ; (I-L) ; (L-I)];
f = [p_l*O ; -p_u*O ; alpha_l*O ; -alpha_u*O ; e ; -e ; Z ; Z];

%% Checking final condition of feasibility

temp = F*p;

for i=1:size(f,1)
    if temp(i)<f(i)
        y=Inf;
        return;
    end
end

%% checking feasibility of slack variable given vector w

for i=1:size(s,1)
    if s(i)>w(i)
        y=Inf;
        return;
    end
end

%% Checking feasibility of equality constraint
temp = zeros(size(w,1),1);
for i=1:N
    temp = temp + Psi{i}*[p(:,i);q(:,i);s(:,i)];
end

for i=1:size(w,1)
    if ~(temp(i)==w(i))
        y=Inf;
        return
    end
end



%% 
% At this point all feasibility conditions should have been met, so we can
% return 0
y=0;

end

