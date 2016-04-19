function [c,d,b_avg,P_modes]=hist_sect(A,P,B=10);
% [c]=hist_sect(A,P,B=10);
% Using data A, find P centres of P regions of equal area of the histogram of A.
% There is a difference if the area is accumulated from the beginning of the
% array or the end. For this reason, the accumulation of both directions is
% computed, the centres of these regions are found and then the centres are
% averaged.
%
% A : The data with which to compute the histogram
% P : The number of centres pursued. Should be less than B.
% B : The number of bins of the histogram, default 10.
%
% c : The centres
% d : The distance of each centre from its closest boundary. This number is
%     always positive.
% b : The boundaries.
% P_modes : Also returns the coordinates of local maxima (the modes).
% extreme bin boundaries
bb_s=min(A);
bb_e=max(A);
% bin boundaries
bbs=(bb_e-bb_s)*(0:B)/B+bb_s;
% bin centres
bcs=(bb_e-bb_s)*((0:(B-1))+0.5)/B+bb_s;
% gather data into bins
[nn,xx]=hist(A,bcs);
[nns,nnsi]=sort(nn,'descend');
P_modes=xx(nnsi(1:P));
P_modes=P_modes(:);
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
c=b_m_avg;
% compute average boundaries
b_avg=(b_+flipud(b_2))/2;
d=min(abs(c-b_avg'),[],2);
b=b_avg;
