function [e]=f2perb(f)
% F2PERB
%
% Convert frequency to perceptual ERB.
%
% Arguments:
%   f: frequency in Hz
e=9.26*log(4.37e-3*f+1);
