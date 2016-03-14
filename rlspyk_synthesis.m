function [W,d] = rlspyk_synthesis(E,p,lambda)
% RLSPYK_SYNTHESIS - RLS predictor in synthesis configuration as implemented by
% Yu and Ko in "Lossless Compression of Digital Audio Using Cascaded RLS-LMS
% Prediction," IEEE Trans. Speech Audio Processing, Vol. 11, No. 6, November
% 2003.
M = length(E);
d = zeros(p+M,1);
delta=0.001;
if nargin < 3, lambda = 1.0; end;
P=eye(p)/delta;
W(1,:)=zeros(1,p);
for k=1:M,
    x=d(p-[-(k-1):-(k-p)]);
    d(p+k)=E(k)+W(k,:)*x;
    z=(1/lambda*P*x)/(1+1/lambda*x.'*P*x);
    W(k+1,:)=W(k,:)+z.'*E(k);
    P=(P-z*x.'*P)/lambda;
end;
d = d((p+1):end); % remove 0 padding
