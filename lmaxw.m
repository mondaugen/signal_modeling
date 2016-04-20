function [ma_i]=lmaxw(x,b,o)
% Find local maxima of function by only looking in a region of size b for a
% maximum.
% 
% x is signal to find maxima in
% b is window length
% o is overlap
%
% ma_i are the indices in x of the maxima.
N=length(x);
ci=(0:o:(N-b));
ri=(1:b)';
S=ri+ci;
[m,mi]=max(x(S));
ma_i=unique(mi+ci);
