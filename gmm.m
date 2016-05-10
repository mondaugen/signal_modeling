function [C,P,mu,S,W]=gmm(X,mu,S,W,B)
%[C,P,mu,S,W]=gmm(X,mu,S,W)
% Gaussian mixture model
% Starting with initial guesses for means mu, covariances S and weights W, use
% the expectation maximization (EM) algorithm to update these parameters and
% also give a classification probability P for each datum in X.
% mu is a matrix of KxN entries where N is the length of each row of X and K is
% the number of classifications.
%
% K is determined from the dimensions of mu.
% S is a tensor of NxNxK entries.
% W is a vector of 1xK entries
% B is the number of iterations of the EM algorithm to compute.
% 
% C are the maximum likelihood classifications
% P are the classification probabilities
% mu are the final means
% S are the final variances
% W are the final weights
K=size(mu,1);
N=size(X,2);
P=zeros(size(X,1),K);
for b=1:B
    for k=1:K
        % In the case the inverse of the variance is Inf, we do not want
        % argument to equal 0 or the function becomes undefined.
        X_mu=X-mu(k,:);
        % Add offset to values equal to 0.
        X_mu(find(X_mu==0))=realmin();
        S_=S(:,:,k);
        S_(find(S_==0))=realmin();
        P(:,k)=W(k)*((2*pi)^N*det(S_))^(-0.5)...
            *exp(-0.5*sum(X_mu'.*(inv(S_)*X_mu'),1)');
        % Add offset to values equal to 0.
        pk0_=find(P(:,k)==0);
        P(pk0_,k)=realmin();
    end
    P=P./sum(P,2);
    W=sum(P,1)/size(P,1);
    for k=1:K
        mu(k,:)=sum(X.*P(:,k),1)/sum(P(:,k),1);
        S(:,:,k)=((X-mu(k,:)).*P(:,k))'*(X-mu(k,:))/sum(P(:,k),1);
    end
end
[c_,C]=max(P,[],2);
