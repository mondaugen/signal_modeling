function [f]=perb2f(e)
% PERB2F
%
% Convert perceptual ERB index to frequency
%
% Arguments:
%   e:  ERB index
f=(exp(e/9.26)-1)/4.37e-3;
