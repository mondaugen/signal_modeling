function [A,E] = lms_lp_analysis(d,mu,p,a0)
% LMS_LP_ANALYSIS -- produces the A coefficients of an adaptive linear predictor
% of order p and its error signal E
% Based on Hayes's implementation of LMS (Hayes p. 506)
% pad the beginning of d with p zeros
d = [zeros(p,1); d(:)];
% make a way to look up negative values in d all the way to -p + 1
d_ = @(n) d((n) + p);
if nargin < 5, a0 = zeros(1,p); end;
a0 = a0(:).';
E(1) = d_(1) - a0*d_(-[0:-(1-p)]);
A(1,:) = a0 + mu*E(1)*conj(d_(-[0:-(1-p)])).';
M = length(d);
if M>1,
    for k=2:M-p,
        E(k) = d_(k) - A(k-1,:)*d_(-[-(k-1):-(k-p)]);
        A(k,:) = A(k-1,:) + mu*E(k)*conj(d_(-[-(k-1):-(k-p)])).';
    end;
end;
