function [lmx,lmxi]=lmax(x)
% [lmx,lmxi]=lmax(x)
% output local maxima and their indices
d=diff(x);
% only consider maxima greater than threshold
d_max=([0;(d>0)]&[(d<0);0]);
lmxi=find(d_max);
lmx=x(lmxi);
