% random sample
X=rand(10000,1);
N_b=500; % number of bins
% bin boundaries
bbs=(max(X)-min(X))*(0:N_b)/N_b+min(X);
bcs=(max(X)-min(X))*((0:(N_b-1))+0.5)/N_b+min(X);
[nn,xx]=hist(X,bcs);
xxi_small=find(nn<2);
bbssu=bbs(1+xxi_small);
bbssl=bbs(xxi_small);
Xdi=find(sum((X<=bbssu)&(X>=bbssl),2));
Y=X;
Y(Xdi)=[];
newplot(figure(1));
newplot(figure(2));
figure(1);
bar(xx,nn);
figure(2);
hist(Y,bcs);
