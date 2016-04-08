function [p_scores]=p_f0_score(f0,p,Ap,max_f,sd_f,h_slp)
% [p_scores]=p_f0_score(f0,p,max_freq,sd_f,harm_dcy)
%
% f0 : fundamental frequency to score
% p : partials to test in scoring
% Ap : the amplitudes of these partials
% max_f : maximum centre frequency of comb tooth
% sd_f : 60 % of the total energy of a harmonic comb tooth is between 
%        f0*k - sd_f and f0*k + sd_f. This may change to increase sd_f with
%        frequency as well.
% h_slp : The rate of change in the amplitudes of the harmonic comb in
%         db/octave. Negative values give a decay.
p=p(:)';
Ap=Ap(:)';
f=f0:f0:max_f;
f=f(:);
dcy=(f./f(1)-1)*h_slp;
Adcy=10.^(dcy/20);
p_scores=(Adcy*Ap).*exp(-((f-p)/sd_f).^2);
p_scores=sum(p_scores,1);
