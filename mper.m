function [Px,ENBW] = mper(x,win,N_fft,L)
% [Px,ENBW] = mper(x,win,N_fft,L)
% 
% Modified periodogram.
%
% Px : The powerspectrum
% ENBW : is the equivalent noise bandwidth summed over all bins of the DFT of
%        the window. This is useful for calculating the total power using Px, as
%        it will have been scaled by this value.
% x : The signal to analyse
% win : can be one of
%       'hamming'
%       'hanning' (default)
%       'bartlett'
%       'blackman'
%       or can be vector of length L, in which case it will be used as the
%       window
% N_fft : length of FFT to perform (default: length of x)
% L : length of window (default: length of x)
%     Window padded with 0s to equal N.
x=x(:);
N=length(x);
if (nargin < 4)
    L=N;
end
L_odd=mod(L,2);
if (nargin == 1)
    win='hanning';
end
if (ischar(win))
    if (strcmp('hanning',win) == 1),
        if L_odd == 1,
            w=hanning(L);
        else
            w=[hanning(L-1); 0];
        end;
    elseif (strcmp(win,'hamming') == 1),
        if L_odd == 1,
            w=hamming(L);
        else
            w=[hamming(L-1); 0];
        end;
    elseif (strcmp(win,'bartlett') == 1),
        if L_odd == 1,
            w=bartlett(L);
        else
            w=[bartlett(L-1); 0];
        end;
    elseif (strcmp(win,'blackman') == 1),
        if L_odd == 1,
            w=blackman(L);
        else
            w=[blackman(L-1); 0];
        end;
    elseif (strcmp(win,'boxcar') == 1),
        w=ones(L,1);
    end;
else
    if (length(win)==L)
        w=win(:);
    else
        error('Bad window length.');
    end
end
if nargin < 3,
    N_fft=N;
end;
if (L > N_fft)
    error('Window length cannot be greater than FFT length.');
end
if (L < N)
    w=[w;zeros(N-L,1)];
end
Px=(abs(fft(x.*w,N_fft))/sum(w)).^2;
ENBW=sum(w.^2)/(sum(w).^2)*N_fft;
