function [R]=acorrm(x,p)
% ACORRM
% Estimate autocorrelation matrix using averaged periodograms
% R is pxp matrix
r=acorr(x,p);
R=toeplitz(r);
