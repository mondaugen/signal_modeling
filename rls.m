function [W,E] = rls(x,d,nord,lambda)
%
delta=0.001;
X=convm(x,nord);
[M,N] = size(X);
if nargin < 4, lambda = 1.0; end;
P=eye(N)/delta;
W(1,:)=zeros(1,N);
for k=2:M-nord+1,
    z=P*X(k,:)';
    g=z/(lambda+X(k,:)*z);
    alpha=d(k)-X(k,:)*W(k-1,:).'; % a priori error
    W(k,:)=W(k-1,:)+alpha*g.';
    E(k)=d(k)-X(k,:)*W(k,:).';    % a posteriori error
    P=(P-g*z.')/lambda;
end;
    
