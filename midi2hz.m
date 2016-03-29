function [f0]=midi2hz(pchs)
% [f0]=midi2hz(pchs)
f0=440*2.^((pchs-69)/12);
