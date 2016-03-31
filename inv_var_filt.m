function [y]=inv_var_filt(x,K,alpha=1)
% [y]=inv_var_filt(x,K,alpha)
N=length(x);
y=zeros(N,1);
for n=1:K
    y(n)=x(n);
end
for n=(K+1):(N-K)
    v=var(x((n-K):(n+K)))*alpha;
    y(n)=x((n-K):(n+K))(:)'*exp(-((-K:K)'/v).^2);
    y(n)/=norm(x((n-K):(n+K))(:))*norm(exp(-((-K:K)'/v).^2));
end
for n=(N-K+1):N
    y(n)=x(n);
end
