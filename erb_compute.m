function [erbs,centres]=erb_compute(ol=0.5,n_max=20)
% ERB_COMPUTE
%
% Compute ERBs.
% The bandwidths go in erbs, the centre frequencies in centres.
% bandwidth and centres are in Hz.
%
% ol is the amount the centre frequency increases as a fraction of the last
% computed bandwith. 1 means the next centre frequecy will go at the edge of the
% last computed bandwidth and less than 1 means it will be centred more closely.
% Typically 0 < ol <= 1.
%
% n_max is the number of bands to compute.
%
% The coefficients to compute the ERBs are from:
% J. Acoust. Soc. Am. 74 (3), September 1983 (Moore and Glasberg)
erbs=zeros(n_max,1);
centres=zeros(n_max,1);
f=0;
a=6.23;
b=93.39;
c=28.52;
for n=1:n_max
    erbs(n)=a*(f/1000.)^2+b*(f/1000.)+c;
    centres(n)=f;
    f+=erbs(n)*ol;
end
