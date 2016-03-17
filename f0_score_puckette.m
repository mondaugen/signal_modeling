function [s]=f0_score_puckette(f0,f,Pf,K=20,sgm=1,test=0,d=@d_default,w=@w_default)
% [s]=f0_score_puckette(f0,f,Pf,K=20,sgm=1,test=1,d=@d_default,w=@w_default)
%
% Determine the "score" for a given f0 for a set of frequecies.
%
% Input arguments:
% f0 : the fundamental to obtain a score for.
% f  : the set of frequencies to use in calculating the score.
% Pf : the power at these frequencies
% K  : number of harmonics to consider
% sgm: paramter for weighting function
% d  : a function that weights frequency scores based on their frequency. Will
%      be called with the arguments k,f0,f_ where k is the harmonic number
%      (starting at 1), f0 is as above, and f_ is a frequency in f. By default
%      d(k,f0,f_)=1/k. f and k can be vectors and in this case f must be a
%      column vector and k must be a row vector.
% w  : a function that is convolved with the weighted harmonic "comb" (weighted
%      by the function d above) that is called with the arguments k,f0,f_,sgm.
%      By default the function is w(k,f0,f_)=1/sqrt(2*pi)*exp((f-k*f0)^2/2). f
%      and k can be vectors and in this case f must be a column vector and k
%      must be a row vector.
if (length(f) == 0)
    s=0;
else
    k=1:K;
    s=w(k,f0,f,sgm)*d(k,f0,f)'.*Pf;
    if (test==0)
        s=sum(s);
    end
end

function [dw]=d_default(k,f0,f)
k=k(:)';
dw=k.^(-1);

function [ww]=w_default(k,f0,f,sgm)
%ww=exp(-(f-k*f0).^(2)/(2*sgm^2))/(sgm*(2*pi)^(1/2));
f=f(:);
k=k(:)';
ww=exp(-(f-k*f0).^(2)/(2*sgm^2));
