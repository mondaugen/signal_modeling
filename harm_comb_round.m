function [c]=harm_comb_round(f,f0)
%function [c]=harm_comb_round(f,f0)
f=f(:);
f0s=f0:f0:max(f);
f0s=f0s(:)';
[mi,ns]=min(abs(f-f0s));
c=zeros(length(f),1);
c(ns)=1;
