function [ex]=mark_extrem(x)
% [ex]=mark_extrem(x)
% remove vector that marks each data point in x as minimum (=-1), maximum (=1),
% or neither (=0).
ex=zeros(length(x),1);
d=diff(x);
d_max=([0;(d>0)]&[(d<0);0]);
d_min=([0;(d<0)]&[(d>0);0]);
ex(find(d_max==1))=1;
ex(find(d_min==1))=-1;
if (x(1)>x(2))
    ex(1)=1;
elseif (x(1)<x(2))
    ex(1)=-1;
end
if (x(end-1)>x(end))
    ex(end)=-1;
elseif (x(end-1)<x(end))
    ex(end)=1;
end
