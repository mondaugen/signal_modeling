function [x]=sum_cos_win_t(n,a,N)
% [x]=sum_cos_win_f(n,a)
% n are the samples at which to evaluate
%   n cannot be a matrix
% a are the coefficients of each cosine in the sum
% N is the length of the DFT that transformed the window (if n are the bin
%   frequencies of an FFT, this is just the length of n)
n=n(:);
if ischar(a)
    a=rcwn2a(a);
end
a=a(:)';
M=length(a);
m=(0:(M-1));
a=a.*((-1).^m);
x=a.*cos(2*pi/N*n*m);
x=sum(x,2);
