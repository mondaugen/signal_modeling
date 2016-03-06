function [a,err]=acm(x,p)
%
x=x(:);
N=length(x);
if p>=length(x)
    error('Model order too large.');
end
X=convm(x,p+1);
Xq=X(1:N+p-1,1:p);
a=[1;linsolve(-Xq,X(2:N+p,1))];
err=abs(X(1:N+p,1)'*X*a);
