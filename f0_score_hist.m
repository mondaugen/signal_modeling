function [s]=f0_score_hist(pchs,f,Pf,alpha,beta)
% [s]=f0_score_hist(f0,f,Pf)
%
% Determine the scores for given pchs compared to a set f of frequencies
%
% Inbut arguments:
% pchs : the pitches to obtain scores for.
%        pitches are such that the f0 in Hz is 440*2^((pch-69)/12)
% f  : the set of frequencies to use in calculating the scores.
% Pf : the power at these frequencies
% alpha : The relative weight of pitches at extreme ends of the pitch table,
%         E.g., alpha = 0.1 will weight pitches equal to min(pchs) and max(pchs)
%         with a weight of 0.1 versus a weight of 1 for pitches equal to
%         (max(pchs)-min(pchs))/2
% beta : The relative weight of partials at a rate per 6dB. E.g, if beta = 0.9
%        and the ratio of partial_i and partial_k's amplitudes are +-6dB, than
%        this comparison will be weighted by 0.9, versus 1 for partials of the
%        same amplitude (ratio of 0 dB).
edges=[pchs(1:(end-1))-diff(pchs)/2 (pchs(2:end)+diff(pchs)/2)((end-1):end)];
s=zeros(length(pchs),1);
if (length(f) > 0)
    fdiffs=f(:)-f(:)';
    Pdiffs=10*log(Pf(:)*((Pf(:)'.^(-1))))/log(10)';
    fdiffs=fdiffs(:);
    Pdiffs=Pdiffs(:);
    pch_min=min(pchs);
    pch_max=max(pchs);
    pch_mid=(pch_min+pch_max)/2;
    pch_diff=pch_max-pch_min;
    sgm_f=(-(pch_diff^2)/(8*log(alpha)))^(1/2);
    sgm_P=(-(6^2)/(2*log(beta)))^(1/2);
    fw=@(x) exp(-(x^2)/(2*sgm_f^2));
    Pw=@(x) exp(-(x^2)/(2*sgm_P^2));
    
    for k=1:length(fdiffs)
        idx=lookup(edges,12*log(abs(fdiffs(k))/440)/log(2)+69);
        if (idx > 0) && (idx <= length(s))
            s(idx)+=fw(fdiffs(k))*Pw(Pdiffs(k));
        end
    end
end
