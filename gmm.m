function [C,P,mu,S,W]=gmm(X,mu,S,W,B)
%[C,P,mu,S,W]=gmm(X,mu,S,W)
% Gaussian mixture model
% Starting with initial guesses for means mu, covariances S and weights W, use
% the expectation maximization (EM) algorithm to update these parameters and
% also give a classification probability P for each datum in X.
% mu is a matrix of KxN entries where N is the length of each row of X and K is
% the number of classifications.
% K is determined from the dimensions of mu.
% S is a tensor of NxNxK entries.
% W is a vector of 1xK entries
% B is the number of iterations of the EM algorithm to compute.
K=size(mu,2);
N=size(X,2);
P=zeros(size(X,1),K);
for b=1:B
    for k=1:K
        P(:,k)=W(k)*((2*pi)^N*det(S(:,:,k)))^(-0.5)...
            *exp(-0.5*sum((x-mu(k,:))'.*(inv(S(:,:,k))*(x-mu(k,:))'),1)');
    end
    P=P./sum(P,2);
    W=sum(P,1)/size(P,1);
    for k=1:K
        mu(k,:)=sum(X.*P(:,k),1)/sum(P(:,k),1);
        S(:,:,k)=((x-mu(k,:)).*P(:,k))'*(x-mu(k,:))/sum(P(:,k),1);
    end
end
C=max(P,[],2);
