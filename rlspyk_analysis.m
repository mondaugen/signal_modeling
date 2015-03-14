function [W,E] = rlspyk_analysis(d,p,lambda)
% RLSPYK_ANALYSIS - RLS predictor in analysis configuration as implemented 
% by Yu and Ko in "Lossless
% Compression of Digital Audio Using Cascaded RLS-LMS Prediction," IEEE
% Trans. Speech Audio Processing, Vol. 11, No. 6, November 2003.
d = [zeros(p,1); d(:)];
delta=0.001;
if nargin < 3, lambda = 1.0; end;
P=eye(p)/delta;
W(1,:)=zeros(1,p);
M = length(d);
for k=1:M-p,
    x=d(p-[-(k-1):-(k-p)]);
    E(k)=d(p+k)-W(k,:)*x;
    z=(1/lambda*P*x)/(1+1/lambda*x.'*P*x);
    W(k+1,:)=W(k,:)+z.'*E(k);
    P=(P-z*x.'*P)/lambda;
end;