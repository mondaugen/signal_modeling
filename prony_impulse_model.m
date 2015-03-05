function [y,err] = prony_impulse_model(x,p,q)
% PRONY_IMPULSE_MODEL -- Models a signal x with p poles and q zeros and outputs
% a signal made by running an impulse thorugh a filter with the estimated
% parameters
N = length(x);
[a,b,err] = prony(x,p,q);
d = zeros(1,N);
d(1)=1;
y = filter(b,a,d);
