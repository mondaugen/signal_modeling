function [W,E,apri,apost] = rls(x,d,nord,lambda)
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
    apri(k)=X(k,:)*W(k-1,:).'
    alpha=d(k)-apri; % a priori error
    W(k,:)=W(k-1,:)+alpha*g.';
    apost(k)=X(k,:)*W(k,:).';
    E(k)=d(k)-apost;    % a posteriori error
    P=(P-g*z.')/lambda;
end;
    
