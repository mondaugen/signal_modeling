function [X]=sum_cos_win_f(theta,a,N)
% [X]=sum_cos_win_f(theta,a)
% theta are frequencies in radians at which to evaluate
% a are the coefficients of each cosine in the sum
% N is the length of the DFT that transformed the window (if theta are the bin
%   frequencies of an FFT, this is just the length of theta)
theta=theta(:);
a=a(:)';
M=length(a);
m=(0:(M-1));
a=a.*((-1).^m);
ofst=2*pi/N*m;
X=0.5*a.*(dk(theta-ofst,N)+dk(theta+ofst,N));
X=sum(X,2);
