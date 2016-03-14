function [Px,ENBW] = mper(x,win,N_fft)
% MPER = Modified periodogram.
%   win can be one of
%    'hamming'
%    'hanning'
%    'bartlett'
%    'blackman'
% If N_fft is given, this will be the length of the FFT that is performed,
% otherwise the length is the length of x.
% ENBW is the equivalent noise bandwidth summed over all bins of the DFT of the
% window. This is useful for calculating the total power using Px, as it will
% have been scaled by this value.
x=x(:);
N=length(x);
N_odd=mod(N,2);
if nargin == 1,
    w=hanning(N);
else
    if (strcmp('hanning',win) == 1),
        if N_odd == 1,
            w=hanning(N);
        else
            w=[hanning(N-1); 0];
        end;
    elseif (strcmp(win,'hamming') == 1),
        if N_odd == 1,
            w=hamming(N);
        else
            w=[hamming(N-1); 0];
        end;
    elseif (strcmp(win,'bartlett') == 1),
        if N_odd == 1,
            w=bartlett(N);
        else
            w=[bartlett(N-1); 0];
        end;
    elseif (strcmp(win,'blackman') == 1),
        if N_odd == 1,
            w=blackman(N);
        else
            w=[blackman(N-1); 0];
        end;
    elseif (strcmp(win,'boxcar') == 1),
        w=ones(N,1);
    end;
end;
if nargin < 3,
    N_fft=N;
end;
Px=(abs(fft(x.*w,N_fft))/sum(w)).^2;
ENBW=sum(w.^2)/(sum(w).^2)*N_fft;
