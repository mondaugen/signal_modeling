function [A,d] = nlms_lp_synthesis(E,be,p,a0)
% NLMS_LP_SYNTHESIS -- produces the A coefficients of an adaptive linear predictor
% of order p and its original signal d from the error signal E 
% Based on Hayes's implementation of LMS (Hayes p. 506)
% pad the beginning of d with p zeros and initialize to 0
M = length(E);
d = [zeros(p+M,1)];
if nargin < 5, a0 = zeros(1,p); end;
a0 = a0(:).';
d(1+p) = E(1) + a0*d(p-[0:-(1-p)]);
D = d(p-[0:-(1-p)])' * d(p-[0:-(1-p)]) + 0.0001;
A(1,:) = a0 + be/D*E(1)*conj(d(p-[0:-(1-p)])).';
if M>1,
    for k=2:M,
        d(k+p) = E(k) + A(k-1,:)*d(p-[-(k-1):-(k-p)]);
        D = d(p-[-(k-1):-(k-p)])' * d(p-[-(k-1):-(k-p)]) + 0.0001;
        A(k,:) = A(k-1,:) + be/D*E(k)*conj(d(p-[-(k-1):-(k-p)])).';
    end;
end;
d = d((p+1):end); % remove 0 padding
