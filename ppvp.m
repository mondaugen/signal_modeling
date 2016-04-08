function [ex_i,X]=ppvp(Y,a=0.1,mode='max',c=0.5,q=2)
% function [ex_i,X]=ppvp(Y,mode='max')
%
% Find set of extrema based on minimizing variation.
%
% Y is data record
% a is recuction factor, e.g., 0.1 means a total variation reduction of 10 % of
%   the original variation will halt the algorithm
% mode can be 'min' or 'max'
%   'min' finds set of minima
%   'max' finds set of maxima
% c is percentage of variation score that is attributed to a values total
%   difference from the maximum, resp. minimum
% q is exponent to which distance from highest resp. lowest value is raised if
%   algorithm is seeking minima resp. maxima. q>1 means two distant values
%   are more different than two close ones and q<1 means the opposite
%
% ex_i are indices in Y that have reached variation reduction threshold
% X are the values of these indices
N=length(Y);
X=Y;
ex_states=mark_extrem(X);
ex_i=(1:N);
va=abs(X-[X(1);X(1:(N-1))])+abs(X-[X(2:N);X(N)]);
tot_va=sum(va);
thresh=tot_va*a;
switch mode
case 'max'
    b=-1;
case 'min'
    b=1;
otherwise
    error(sprintf('Bad mode %s',mode));
end
while ((tot_va > thresh) && (N>1))
    mi_i=find(ex_states==b);
    if (length(mi_i)==0)
        % No more items with state b exist
        break
    end
    [va_ma,va_ma_i]=max(va(mi_i));
    va_ma_i=ex_i(mi_i(va_ma_i));
    ex_i=ex_i(find(ex_i != va_ma_i));
    X=Y(ex_i);
    ex_states=mark_extrem(X);
    N-=1;
    switch mode
    case 'max'
        X_ex=max(X);
    case 'min'
        X_ex=min(X);
    end
    va=(1-c)*abs(X-[X(1);X(1:(N-1))])+abs(X-[X(2:N);X(N)]);
    tot_va=sum(va);
    va+=c*abs(X-X_ex).^(q);
end
