function [W,E] = rls_lp_analysis(d,p,lambda)
% RLS_LP_ANALYSIS -- produces the A coefficients of an adaptive linear predictor
% of order p and its error signal E using the Recursive Least-Squares algorithm
% Based on Hayes's implementation of RLS (Hayes p. 546)
% pad the beginning of d with p zeros
d = [zeros(p,1); d(:)];
delta=0.001;
if nargin < 3, lambda = 1.0; end;
P=eye(p)/delta;
W(1,:)=ones(1,p);
M = length(d);
for k=2:M-p,
    z=P*d(p-[-(k-1):-(k-p)]);
    g=z/(lambda+d(p-[-(k-1):-(k-p)]).'*z);
    alpha=d(k+p)-W(k-1,:)*d(p-[-(k-1):-(k-p)]); % a priori error
    W(k,:)=W(k-1,:)+alpha*g.';
    E(k)=d(k+p)-W(k,:)*d(p-[-(k-1):-(k-p)]); % a priori error
    P=(P-g*z.')/lambda;
end;
