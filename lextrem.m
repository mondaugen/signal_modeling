function [lmx,lmxi]=lextrem(x,mode='max')
% [lmx,lmxi]=lextrem(x,dir='max')
% output local extrema and their indices
% mode can be 'min' or 'max'
d=diff(x);
switch mode
case 'max'
    d_max=([0;(d>0)]&[(d<0);0]);
case 'min'
    d_max=([0;(d<0)]&[(d>0);0]);
otherwise
    error(sprintf('Bad mode %s.',mode));
end
lmxi=find(d_max);
lmx=x(lmxi);
