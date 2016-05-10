function [y]=ac_filt(x,b,dim=1)
% [Y]=AC_FILT(X,B)
%
% Filter X with the anticausal filter B. That means (with 0 based indexing you
% shitheads)
%   y[n] = sum_(k=-floor(length(b)/2)^(floor(length(b)/2))
%                       x[n-k]*b[floor(len(b)/2)+k]
% B must have an odd length.
% Filters x along dimension "dim" (default 1)
N=size(b);
if (mod(b,2)==0)
    error(sprintf('B must have odd length. Length is %d',N))
end
N_2=floor(N/2);
switch (dim)
case 1
    x=[x;zeros(N_2,size(x,2))];
case 2
    x=[x,zeros(size(x,1),N_2)];
end
y=filter(b,1,x,[],dim);
switch (dim)
case 1
    y=y((N_2+1):end,:);
case 2
    y=y(:,(N_2+1):end);
end
