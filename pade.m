function [a,b] = pade(x,p,q)
% PADE -- see Hayes p. 138
x = x(:);
if p+q>=length(x), error('Model order too large.'), end;
X = convm(x,p+1);
Xq = X(q+2:q+p+1,2:p+1);
a = [1;-Xq\X(q+2:q+p+1,1)];
b = X(1:q+1,1:p+1)*a;
end;
