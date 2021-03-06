function [a,b,err] = prony(x,p,q)
% PRONY -- See Hayes p. 154
x = x(:);
N = length(x);
if p+q>=length(x),
    error(strcat('Model order too large. ',...
        ' length= ',num2str(length(x)),...
        ' p= ',num2str(p),...
        ' nq= ',num2str(q)));
end;
X = convm(x,p+1);
Xq = X(q+1:N+p-1,1:p);
a = [1;-Xq\X(q+2:N+p,1)];
%a = [1;linsolve(X(q+2:N+p,1),-Xq)];
b = X(1:q+1,1:p+1)*a;
err = x(q+2:N)'*X(q+2:N,1:p+1)*a;
end;
