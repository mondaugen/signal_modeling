% Realize random variable as mixture of two different distributions
mu1=1;
mu2=10;
sig1=.3;
sig2=.1;
p1=0.5;
% Number of realizations
N=400;
R_sp=0.5;
sp_min=-20;
B=500;
sp_max=20;
D=unifrnd(0,1,N,1);
X=(D<=p1).*(randn(N,1)*sig1+mu1)+(D>p1).*(randn(N,1)*sig2+mu2);
% Add spurious data
X=[X;unifrnd(sp_min,sp_max,round(N*R_sp),1)];
% Shuffle
X=X(randperm(length(X)));
X=X(:);
figure(1);
hist(X,100);
S=zeros(1,1,2);
S(1,1,:)=[1;1];
[C,P,mu,S,W]=gmm(X,[0;10],S,[0.5,0.5],B);
%[C,P,mu,S,W]=gmmr(X,[0;10],S,[0.5,0.5],B,100);
newplot(figure(2));
figure(2);
hold on;
scatter(X(find(C==1)),zeros(length(find(C==1)),1),[],'r');
scatter(X(find(C==2)),zeros(length(find(C==2)),1),[],'g');
hold off;
