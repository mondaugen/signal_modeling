% Find the P centres of parts of equal area of a histogram.
% This doesn't work but it is the right idea. It doesn't work because interp1
% will not work with multiple equal x values.
P=7;
A_mi=0.9;
A_ma=15;
N=100;
%A=unifrnd(A_mi,A_ma,N,1);
A=A_mi+randn(2*N,1);
A=[A;A_ma+randn(N,1)];
% extreme bin boundaries
bb_s=min(A);
bb_e=max(A);
% number of bins
B=10;
% bin boundaries
bbs=(bb_e-bb_s)*(0:B)/B+bb_s;
% bin centres
bcs=(bb_e-bb_s)*(1:B)/(B+1)+bb_s;
% gather data into bins
[nn,xx]=hist(A,bcs);
cs=cumsum(nn);
i_=interp1([0,cs]/cs(end),(0:length(cs))/length(cs),(1:P)/(P+1))*length(cs);
xx_=interp1(0:length(cs),bbs,i_);
hist(A,bcs);
hold on;
stem(xx_,max(nn)*ones(1,length(xx_)),'r');
hold off;
