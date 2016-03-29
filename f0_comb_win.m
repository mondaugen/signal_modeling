function [w1]=f0_comb_win(x,mu,t,k,p,a)
% function [w]=f0_comb_win(x,mu,t,k,p,a)
% A comb-like window with a peak at mu and negative peaks at periods t.
% k controls the sharpness, while p controls the number of periods before the
% overall window becomes 0, which controls the rate of decay of the height of
% the peaks.
% a controls the width of a gaussian window windowing all the values relative to
% the number of periods p
%w=1-abs((x-mu)/(p*t));
%w=w.*(1.5*exp(-(sqrt(k)*(x-mu)/t).^2)+(1+exp(-2*k*(cos(pi*((x-mu)/t-0.5))).^2)).^(-1)-1);
%w=w.*((1+exp(-2*k*(cos(pi*((x-mu)/t-0.5))).^2)).^(-1)-1);
%w=-(exp(-k*(cos(pi*((x-mu)/t-0.5))).^2));
x=x(:);
mu=mu(:)';
t=t(:)';
p=p(:)';
k=k(:)';
a=a(:)';
alpha=-k*16./(2*t.^2);
%w1=(1-abs((x-mu)/(p*t)));
w1=exp(-((x.-mu)./(a.*p.*t/2)).^2).*((2*exp(alpha.*(x.-mu).^2)-exp(-k.*(cos(pi*((x.-mu)./t-0.5))).^2)).+exp(-k))./(1.+exp(-k));
%w2=alpha*(x-mu).^2;
%w3=pi*((x-mu)/t-0.5);
