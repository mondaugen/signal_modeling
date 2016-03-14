function [B]=div_ival_ol(a,N,bet)
% DIV_IVAL_OL
%
% Divide length a into N equally sized intervals overlapping by bet.
%
% Arguments:
% a: length
% N: number of intervals
% bet: overlapping factor (0-1)

c=(0:(a/N):(a-a/N))+a/(2*N);
B=zeros(length(c),2);
B(:,1)=c-a/N*bet;
B(:,2)=c+a/N*bet;
