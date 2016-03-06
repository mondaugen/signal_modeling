function [r]=acorr(x,p)
% ACORR
% Estimate autocorrelation using averaged periodograms
% r is column vector of length p where r(1) = ac_x(0) r(2) = ac_x(1) etc.
x=x(:);
nidx=1:(length(x)-p+1);
r=zeros(p,1);
for n=nidx
    S=fft(x(n:(n+p-1)));
    r+=real(ifft(conj(S).*S/p));
end
r/=length(nidx);
