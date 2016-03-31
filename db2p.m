function [p]=db2p(P)
% [p]=db2p(P)
% convert db Power into linear power
p=10.^(P/10);
