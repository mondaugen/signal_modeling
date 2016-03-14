function [a,err]=covm(x,p)
%
x=x(:);
N=length(x);
if p>=length(x)
    error('Model order too large.');
end
X=convm(x,p+1);
Xq=X(p:N-1,1:p);
a=[1;linsolve(-Xq,X(p+1:N,1))];
err=abs(X(p+1:N,1)'*X(p+1:N,:)*a);
