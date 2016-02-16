function [w,P,e] = rls_1(d,x,w,P,lamb)
% RLS_1
% One iteration of the RLS algorithm for adaptively determining filter
% coefficients.
% Adapt coefficients w to minimize error in predicting d from w'*x
% w(k) coefficient corresponding to x(end-k+1)
% K in [1,...,p]
x=flipud(x(:));
w=w(:);
z=P*x;
g=z./(lamb+x'*z);
e=d-w'*x;
w=w+e*g;
P=(P-g*z')./lamb;
