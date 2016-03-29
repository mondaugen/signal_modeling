function [X]=sum_cos_win_f(theta,a,N)
% [X]=sum_cos_win_f(theta,a)
% theta are frequencies in radians at which to evaluate
%   theta cannot be a matrix
% a are the coefficients of each cosine in the sum
% N is the length of the DFT that transformed the window (if theta are the bin
%   frequencies of an FFT, this is just the length of theta)
theta=theta(:);
if ischar(a)
    a=rcwn2a(a);
end
a=a(:)';
M=length(a);
m=(0:(M-1));
a=a.*((-1).^m);
ofst=2*pi/N*m;
%X=0.5*a.*(dk(theta-ofst,N)+dk(theta+ofst,N));
%X=0.5*a.*(exp(-j*(N-1)*(theta-ofst)/2).*dk(theta-ofst,N-1)+exp(-j*(N-1)*(theta+ofst)/2).*dk(theta+ofst,N-1));
% Multiply by complex exponential to represent Fourier transform of time-domain window whose
% centre is at N/2
X=0.5*a.*(exp(-j*(N)*(theta-ofst)/2).*dk(theta-ofst,N)+exp(-j*(N)*(theta+ofst)/2).*dk(theta+ofst,N));
X=sum(X,2);
