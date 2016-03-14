function [W]=perb_hann(f,N,bet)
% PERB_HANN
%
% Evaluate N Hann windows that divide the perceptual ERB spectrum into N equal
% bands that overlap by bet at frequencies f
%
% Arguments:
%   f:      the frequencies at which to evaluate
%   N:      the number of divisions of the spectrum and therefore the number of
%           Hann windows
%   bet:    the overlapping factor
e_max=f2perb(max(f));
B=div_ival_ol(e_max,N,bet);
ec=(B(:,1)+B(:,2))/2;
B=(-B(:,1)+B(:,2));
B=B(1);
W=erb_hann(f2perb(f),ec,B);
