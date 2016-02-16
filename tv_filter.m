function [y]=tv_filter(B,A,x)
% TV_FILTER
% 
% Similar to "filter" but B and A must have same number of rows as length of x
% For each index n of x, B(n,:) and A(n,:) are used in the difference equation.
% Also, can only filter one column of x at a time, if you don't pass in a
% single column or row vector, funny things will happen.
% x is forced to be a column vector.
x=x(:);
if (size(B,1) != length(x)) || (size(A,1) != length(x))
    error('B and A must have same number of rows as length of x.');
end
y=zeros(length(x),1);
[y(1),state]=filter(B(1,:),A(1,:),x(1));
for n=2:length(x)
    [y(n),state]=filter(B(n,:),A(n,:),x(n),state);
end
