function [c]=harm_comb_gauss(f,f0,sgm)
% function [c]=harm_comb_gauss(f,f0,sgm)
% f are frequencies at which to evaluate
% f0 is fundamental and frequency of first harmonic
% sgm is the difference in frequency where the function centred on each harmonic
%   attains a value of 1/e of its maximum
f=f(:);
f0s=f0:f0:max(f);
f0s=f0s(:)';
c=exp(-((f-f0s)/sgm).^2);
c=sum(c,2);
