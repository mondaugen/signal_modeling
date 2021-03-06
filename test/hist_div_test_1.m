% Find the P centres of parts of equal area of a histogram.
clear;
P=2;
A_mi=0.9;
A_ma=15;
N=100;
%A=unifrnd(A_mi,A_ma,N,1);
A=A_mi+randn(N,1);
A=[A;A_ma+randn(N,1)];
% extreme bin boundaries
bb_s=min(A);
bb_e=max(A);
% number of bins
B=10;
% bin boundaries
bbs=(bb_e-bb_s)*(0:B)/B+bb_s;
% bin centres
bcs=(bb_e-bb_s)*((0:(B-1))+0.5)/B+bb_s;
% gather data into bins
[nn,xx]=hist(A,bcs);
% cumulate over bins
cs=[0,cumsum(nn)];
cs/=cs(end);
Fc=(1:(P))/(P);
Fc=Fc(:);
% Find indices that just exceed thresholds
nr=nonzeros((cumsum((cs>=Fc),2)==1).*(1:length(cs)));
cs=cs(:);
nl=nr-1;
n_=(Fc-cs(nl))./(cs(nr)-cs(nl))+nl;
% get data based indices
b_=interp1(1:length(cs),bbs,n_);
% add first bin boundary
b_=[bbs(1);b_];
b_m=b_(1:(end-1))+diff(b_)/2;
%% now go backwards to bias the other way
% cumulate over bins
nn2=fliplr(nn);
bbs2=fliplr(bbs);
cs2=[0,cumsum(nn2)];
cs2/=cs2(end);
Fc2=(1:(P))/(P);
Fc2=Fc2(:);
% Find indices that just exceed thresholds
nr2=nonzeros((cumsum((cs2>=Fc2),2)==1).*(1:length(cs2)));
cs2=cs2(:);
nl2=nr2-1;
n_2=(Fc2-cs2(nl2))./(cs2(nr2)-cs2(nl2))+nl2;
% get data based indices
b_2=interp1(1:length(cs2),bbs2,n_2);
% add first bin boundary
b_2=[bbs2(1);b_2];
b_m2=b_2(1:(end-1))+diff(b_2)/2;
% compute average centres
b_m_avg=(b_m+flipud(b_m2))/2;
newplot(figure(1));
hold on
bar(xx,nn);
plot(bbs,ones(length(bbs),1),'or',bcs,ones(length(bcs),1),'ob');
stem(b_(:),max(nn)*ones(length(b_),1),'color','k');
stem(b_m(:),max(nn)*ones(length(b_m),1),'color','r');
stem(b_2(:),max(nn)*ones(length(b_2),1),'color','cyan','marker','x');
stem(b_m2(:),max(nn)*ones(length(b_m2),1),'color','magenta','marker','x');
stem(b_m_avg(:),max(nn)*ones(length(b_m_avg),1),'color','blue','marker','+');
hold off;
