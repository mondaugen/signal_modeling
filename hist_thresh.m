function [Y,Xdi,Xki,bcs,bbs]=hist_thresh(X,B,th)
% function [Y,Xdi,Xki,bcs,bbs]=hist_thresh(X,B,th)
% 
% X : Random sample, column vector
% B : Number of bins
% th : Threshold
%
% Y : Samples that fell into bins containing a number of values greater or equal
%     to th
% Xdi : Indices of bins to be discarded
% Xki : Indices of bins to be kept
% bcs : bin centres
% bbs : bin boundaries

% bin boundaries
bbs=(max(X)-min(X))*(0:B)/B+min(X);
bcs=(max(X)-min(X))*((0:(B-1))+0.5)/B+min(X);
% Guarantee the endpoints are exactly data values
bbs(1)=min(X);
bbs(end)=max(X);
[nn,xx]=hist(X,bcs);
xxi_small=find(nn<th);
bbssu=bbs(1+xxi_small);
bbssl=bbs(xxi_small);
Xdi_v=sum((X<=bbssu)&(X>=bbssl),2);
Xdi=find(Xdi_v);
Xki=find(!Xdi_v);
Y=X;
Y(Xdi)=[];
