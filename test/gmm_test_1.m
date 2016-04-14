% Compute the Principle Components of a dataset using principle components
% analysis (PCA) and then use try and use Gaussian Mixture Models to
% classify the data sets.
clear;
fi=fopen('data/irises.dat','r');
dat=fscanf(fi,'%f ',Inf);
fclose(fi);
dat_size=[12 50];
dat=reshape(dat,dat_size);
dat=dat';
I1=dat(:,1:4);
I2=dat(:,5:8);
I3=dat(:,9:12);
I=[I1;I2;I3];
% The true classification of each datum
c=[ 1*ones(dat_size(1),1);
    2*ones(dat_size(1),1);
    3*ones(dat_size(1),1); ];
% compute principle components
%S=cov(I);
S=corr(I);
[V,L]=eig(S);
[l,li]=sort(diag(L),'descend');
V=V(:,li);
% the PCs
Y1=I1*V;
Y2=I2*V;
Y3=I3*V;
Y=I*V;
newplot(figure(1));
figure(1);
hold on;
h=scatter(Y1(:,1),Y1(:,2),[],'r');
set(h,'linewidth',1);
h=scatter(Y2(:,1),Y2(:,2),[],'g');
set(h,'linewidth',1);
h=scatter(Y3(:,1),Y3(:,2),[],'b');
set(h,'linewidth',1);
hold off;
% The EM algorithm to try and categorize the PCs into the correct groups.
% First we try with the 1st PC only
y=Y(:,1);
% Initial guesses
mu=[3,6,8];
sig=[1,1,1];
W=[1,1,1]/3;
B=2;
for b=1:B
    % Compute the probability y_n is from distribution k
    P_k=(sig*sqrt(2*pi)).^(-1).*exp(-((y-mu).^2)./(2*sig.^2)).*W;
    P_k=P_k./sum(P_k,2);
    % Recompute weights
    mu_old=mu;
    sig_old=sig;
    W_old=W;
    W=sum(P_k,1)/size(P_k,1);
    mu=sum(y.*P_k,1)./sum(P_k,1);
    sig=sum(((y-mu).^2).*P_k,1)./sum(P_k,1);
    hlt=sum(abs(1 - [mu sig W]./[mu_old sig_old W_old]))...
            /(length(W)+length(mu)+length(sig));
    if (hlt < 0.01)
        break;
    end
end
[m,mi]=max(P_k,[],2);
newplot(figure(2));
figure(2)
hold on;
h=scatter(Y(find(mi==1),1),Y(find(mi==1),2),[],'r');
set(h,'linewidth',1);
h=scatter(Y(find(mi==2),1),Y(find(mi==2),2),[],'g');
set(h,'linewidth',1);
h=scatter(Y(find(mi==3),1),Y(find(mi==3),2),[],'b');
set(h,'linewidth',1);
hold off;
% Now we try with the first 2 components
y2=Y(:,1:2);
mu2=[
    3 6 8;
    5.5 5 5.75;
];
sig2=zeros(2,2,3);
sig2(:,:,1)=eye(2);
sig2(:,:,2)=eye(2);
sig2(:,:,3)=eye(2);
W2=[1 1 1]/3;
B=5;
P_k2=zeros(size(y2,1),3);
for b=1:B
    for _k=1:3
        P_k2(:,_k)=((2*pi)^2*det(sig2(:,:,_k)))^(-1)...
            * exp(-0.5*sum((y2-mu2(:,_k)')'.*((sig2(:,:,_k)^(-1))*(y2-mu2(:,_k)')'),1)');
    end
    P_k2=P_k2./sum(P_k2,2);
    mu2_old=mu2;
    sig2_old=sig2;
    W2_old=W2;
    W2=sum(P_k2,1)/size(P_k2,1);
    for _k=1:3
        mu2(:,_k)=(sum(y2.*P_k2(:,_k),1)./sum(P_k2(:,_k),1))';
        sig2(:,:,_k)=((y2-mu2(:,_k)').*P_k2(:,_k))'*(y2-mu2(:,_k)')./sum(P_k2(:,_k));
    end
end
[m2,mi2]=max(P_k2,[],2);
newplot(figure(3));
figure(3)
hold on;
h=scatter(Y(find(mi2==1),1),Y(find(mi2==1),2),[],'r');
set(h,'linewidth',1);
h=scatter(Y(find(mi2==2),1),Y(find(mi2==2),2),[],'g');
set(h,'linewidth',1);
h=scatter(Y(find(mi2==3),1),Y(find(mi2==3),2),[],'b');
set(h,'linewidth',1);
hold off;
