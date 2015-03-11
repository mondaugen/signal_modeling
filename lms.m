function [A,E] = lms(x,d,mu,nord,a0)
%
X = convm(x,nord);
[M,N] = size(X);
if nargin < 5, a0 = zeros(1,N); end;
a0 = a0(:).';
E(1) =  d(1) - a0*X(1,:).';
A(1,:) = a0 + mu*E(1)*conj(X(1,:));
if M>1,
    for k=2:M-nord+1,
        E(k) = d(k) - A(k-1,:)*X(k,:).';
        A(k,:) = A(k-1,:) + mu*E(k)*conj(X(k,:));
    end;
end;
