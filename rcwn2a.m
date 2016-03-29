function [a]=rcwn2a(wname)
% [a]=rcwn2a(wname)
% Convert the name of a raised cosine window into coefficients for its
% computation
% See Harris (1978) for the coefficients
switch wname
case {'blackman3','blackman'}
    a=[0.42323,0.49755,0.07922];
case {'blackman4'}
    % This gives -74dB sidelobes, error in Harris ?
    a=[.35875,.48829,.14128,.01168];
%    a=[.40217,.49703,.09392,.00183];
case {'hann','hanning'}
    a=[0.5,0.5];
case {'hamming'}
    a=[0.54,0.46];
otherwise
    error(sprintf('Bad window: %s',wname));
end
