function [P]=a2db(a)
% [P]=a2db(a)
% convert amplitude to dB Power
P=20*log(a)/log(10);

