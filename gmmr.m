function [C,P,mu,S,W]=gmmr(X,mu,S,W,B,bet)
%[C,P,mu,S,W]=gmmr(X,mu,S,W)
% Regularized Gaussian mixture model
% Starting with initial guesses for means mu, covariances S and weights W, use
% the expectation maximization (EM) algorithm to update these parameters and
% also give a classification probability P for each datum in X.
% mu is a matrix of KxN entries where N is the length of each row of X and K is
% the number of classifications.
%
% X is a matrix of MxN (M observations on N random variables)
% mu is a matrix of KxN (K guesses for the initial means of X)
% K is determined from the dimensions of mu.
% S is a tensor of NxNxK entries.
% W is a vector of 1xK entries
% B is the number of iterations of the EM algorithm to compute.
% bet is the regularization parameter (cannot be 0, must be positive (?))
% 
% C are the maximum apriori classifications
% P are the classification probabilities
% mu are the final means
% S are the final variances
% W are the final weights
K=size(mu,1);
N=size(X,2);
M=size(X,1);
P=zeros(size(X,1),K);
nu=0;
alpha=(N+1)/2;
gam=ones(K,1);
Bet=bet*eye(N);
% Initial mean guess
mu0=mu;
for b=1:B
    for k=1:K
        % In the case the inverse of the variance is Inf, we do not want
        % argument to equal 0 or the function becomes undefined.
        X_mu=X-mu(k,:);
        % Add offset to values equal to 0.
        X_mu(find(X_mu==0))=1e-10;
        P(:,k)=W(k)*((2*pi)^N*det(S(:,:,k)))^(-0.5)...
            *exp(-0.5*sum(X_mu'.*(inv(S(:,:,k))*X_mu'),1)');
    end
    P=P./sum(P,2);
    W=(sum(P,1)+gam'-1)/(M+sum(gam)-N);
    for k=1:K
        mu(k,:)=(sum(X.*P(:,k),1)+nu*mu0(k,:)')/(sum(P(:,k),1)+nu);
        S(:,:,k)=(((X-mu(k,:)).*P(:,k))'*(X-mu(k,:))...
            +nu*(mu(k,:)-mu0(k,:))'*(mu(k,:)-mu0(k,:))+2*Bet)/(sum(P(:,k),1)+2*alpha-N);
    end
end
[c_,C]=max(P,[],2);
