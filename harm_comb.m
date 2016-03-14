function [C]=harm_comb(f,f0,alpha)
% HARM_COMB
%
% Evaluate harmonic comb at frequencies f and with fundamental or period f0.
% alpha controls the steepness, values near 1 are very sharp and values near 0
% are not sharp
f=f(:);
C=(1./(1-alpha*cos(2*pi*f/f0))-(1+alpha)^-1)/((1-alpha)^(-1)-(1+alpha)^(-1));
