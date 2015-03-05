function [a,b,err] = prony(x,p,q)
% PRONY -- See Hayes p. 154
x = x(:);
N = length(x);
if p+q>=length(x), error('Model order too large.'), end;
X = convm(x,p+1);
Xq = X(q+1:N+p-1,1:p);
a = [1;-Xq\X(q+2:N+p,1)];
b = X(1:q+1,1:p+1)*a;
err = x(q+2:N)'*X(q+2:N,1:p+1)*a;
end;
