function [a,hslp]=rcwn2a(wname)
% [a,hslp]=rcwn2a(wname)
% Convert the name of a raised cosine window into coefficients for its
% computation
% a : coefficients
% hslp : highest side-lobe level in dB power
% See Harris (1978) for the coefficients
switch wname
case {'blackman3-min','blackman'}
    a=[0.42323,0.49755,0.07922];
    hslp=-67;
case {'blackman3'}
    a=[0.44959,0.49364,0.05677];
    hslp=-61;
case {'blackman4-min'}
    a=[.35875,.48829,.14128,.01168];
    hslp=-92;
case {'blackman4'}
    a=[.40217,.49703,.09392,.00183];
    hslp=-74;
case {'hann','hanning'}
    a=[0.5,0.5];
    hslp=-32;
case {'hamming'}
    a=[0.54,0.46];
    hslp=-43;
otherwise
    error(sprintf('Bad window: %s',wname));
end
