function [x,X,z,Z,l,L] = ADMM(iteration,n,p,rho,costs)
% Classis ADMM of Boyd's book. Used for comparison
%   Algorithm is found in boyd's admm book, page 49.


x = zeros(p,n);
l = zeros(p,n);
z = zeros(p,1);

X={};
Z={};
L={};

for k=1:iteration
    for i=1:n
    fun = @(x)costs{i}(x) + l(:,i)'*(x-z) + (rho/2)*(x-z)'*(x-z);
    x(:,i) = fminsearch(fun, zeros(p,1));
    end
    
    X{k} = x;

    temp = 0;
    for i=1:n
        temp = temp+ x(:,i) + (1/rho)*l(:,i);
    end
    z = (1/n)*temp;

    Z{k} = z;
    
    for i=1:n
        l(:,i) = l(:,i) + rho*(x(:,i) - z);
    end
    
    L{k} = l;
    
end





end

