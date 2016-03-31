function [P]=p2db(p)
% [P]=a2db(a)
% convert linear power to dB Power
P=10*log(p)/log(10);

