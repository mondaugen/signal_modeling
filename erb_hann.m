function [W]=erb_hann(f,fc,B)
% ERB_HANN
%
% Evaluate the Hann window with centre frequency fc and ERB B at frequencies f.
% In the case that fc and B are vectors, the frequencies are evaluated along the rows
% and the different centre frequencies fc are along the columns
f=f(:);
fc=fc(:)';
B=B(:)';
W=cos(pi./(2*B).*(f-fc)).^2.*((-B<=(f-fc))&((f-fc)<=B));
