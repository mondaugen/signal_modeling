function [Px,ENBW]=welch(x,L,O,win,N_fft)
% WELCH - Average modified periodograms.
% x is the signal whose power spectrum will be estimated
% L is the length of the sub-windows
% O is the amount of overlap of the sub-windows
% win is the kind of window to use. It can be any of mper's accepted windows.
% N_fft is the size of the FFT to perform on each sub-window. If unspecified, it
%   will be set to N=length(x)
% ENBW is the same as in mper (see mper).
x=x(:);
N=length(x);
if (nargin < 5),
    N_fft=N;
end;
end_idxs=(L:(L-O):N);
Px=zeros(N_fft,1);
for e=end_idxs,
    [Px_,ENBW]=mper(x(e-L+1:e),win,N_fft);
    Px+=Px_;
end;
Px/=length(end_idxs); % ceil((N-L)/(L-O))
