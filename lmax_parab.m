function [x_,n_,p]=lmax_parab(x,n)
% [x_,n_]=lmax_parab(x,n)
% Find extreme points of interpolating parabola of the 2 neighbouring points and
% n of x.
x=x(:).';
n=n(:)';
k=1;
m=[-1,0,1];
p=zeros(length(n),3);
for l=n
    if (l==1)
        p(k,:)=polyfit(m,[x(end),x(l:(l+1))],2);
    elseif (l==length(x))
        p(k,:)=polyfit(m,[x((l-1):l),x(1)],2);
    else
        p(k,:)=polyfit(m,x((l-1):(l+1)),2);
    end
    k+=1;
end
n_=-p(:,2)./(2*p(:,1));
x_=p(:,1).*n_.^2+p(:,2).*n_+p(:,3);
n_=n_+n(:);
