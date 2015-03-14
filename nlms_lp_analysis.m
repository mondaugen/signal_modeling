function [A,E] = nlms_lp_analysis(d,be,p,a0)
% NLMS_LP_ANALYSIS -- produces the A coefficients of an adaptive linear predictor
% of order p and its error signal E
% Based on Hayes's implementation of LMS (Hayes p. 506)
% pad the beginning of d with p zeros
% By Nicholas Esterer
d = [zeros(p,1); d(:)];
% make a way to look up negative values in d all the way to -p + 1
d_ = @(n) d((n) + p);
if nargin < 5, a0 = zeros(1,p); end;
a0 = a0(:).';
E(1) = d_(1) - a0*d_(-[0:-(1-p)]);
D = d(p-[0:-(1-p)])' * d(p-[0:-(1-p)]) + 0.0001;
A(1,:) = a0 + be/D*E(1)*conj(d_(-[0:-(1-p)])).';
M = length(d);
if M>1,
    for k=2:M-p,
        E(k) = d_(k) - A(k-1,:)*d_(-[-(k-1):-(k-p)]);
        D = d(p-[-(k-1):-(k-p)])' * d(p-[-(k-1):-(k-p)]) + 0.0001;
        A(k,:) = A(k-1,:) + be/D*E(k)*conj(d_(-[-(k-1):-(k-p)])).';
    end;
end;
