function [W,d] = rls_lp_synthesis(E,p,lambda)
% RLS_LP_ANALYSIS -- produces the W coefficients of an adaptive linear predictor
% of order p and its original signal d from an error signal E using the
% Recursive Least-Squares algorithm Based on Hayes's implementation of RLS
% (Hayes p. 546)
% pad the beginning of d with p zeros
% THIS DOESN'T WORK
M = length(E);
d = zeros(p+M,1);
delta=0.001;
if nargin < 3, lambda = 1.0; end;
P=eye(p)/delta;
W(1,:)=ones(1,p);
for k=2:M,
    z=P*d(p-[-(k-1):-(k-p)]);
    g=z/(lambda+d(p-[-(k-1):-(k-p)]).'*z);
    alpha=d(k+p)-W(k-1,:)*d(p-[-(k-1):-(k-p)]); % a priori error
    W(k,:)=W(k-1,:)+alpha*g.';
    % synthesize from the a posteriori error signal
    d(k+p)=E(k)+W(k,:)*d(p-[-(k-1):-(k-p)]);
    P=(P-g*z.')/lambda;
end;
d = d((p+1):end); % remove 0 padding
