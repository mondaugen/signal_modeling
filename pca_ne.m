function [A,l]=pca_ne(X,mode='cov')
% [A,l]=pca_ne(X,mode='cov')
% Compute the principle components of A and their variances l.
% If mode is 'cov', the covariance matrix of X is used,
% if mode is 'corr', the correlation matrix of X is used.
% X is (NxP) and a matrix of N realizations on P variables.
switch mode
case 'cov'
    S=cov(X);
case 'corr'
    S=corr(X);
otherwise
    error(sprintf('Bad mode %s\n',mode));
end
[V,L]=eig(S);
[l,li]=sort(diag(L),'descend');
V(:,li);
A=X*V;
