function [S,T,F]=stft_ne(s,N,H,W,L,Fs=44100)
% function [S,T,F]=stft_ne(s,N,H,W,L,Fs=44100)
if (L > N)
    error('Window sizes greater than N not yet supported.');
end
switch (W)
    case 'hanning'
        w=hanning(L+1);
    case 'hamming'
        w=hamming(L+1);
    case 'blackman'
        w=blackman(L+1);
    case 'bartlett'
        w=bartlett(L+1);
    case 'rectangular'
        w=onex(L+1,1);
    otherwise
        error(sprintf('Bad window: %s',W));
end
w=w(1:L);
% pad with zeros
w=[w;zeros(N-L,1)];
W_0=sum(w);
% space to load in samples
x=zeros(N,1);
% space to write spectrogram
S=[];
% space to write time
T=[];
% current time at start of window.
n=1;
% Load in enough samples to fill half of x so that the first frame's centre is
% aligned with the first sample
L_end=length(s);
if (L_end < N/2)
    x((N/2+1):(N/2+L_end))=s;
    x((N/2+1+L_end):N)=zeros(N/2-L_end,1);
else
    x((N/2+1):N)=s(1:N/2);
end
n=n+H;
while (n <= length(s))
    T=[T (n-1)/Fs];
    X=fft(x.*w)/W_0;
    S=[S X];
    % shift over values
    x(1:(N-H))=x((H+1):N);
    L_end=length(s)-n+1;
    if (L_end < H)
        x((N-H+1):(N-H+L_end))=s(n:end);
    else
        x((N-H+1):N)=s(n:(n+H-1));
    end
    n=n+H;
end
% pad with zeros until buffer full
x((N-H+1+L_end):N)=zeros(H-L_end,1);
X=fft(x.*w)/W_0;
S=[S X];
F=(0:(size(S,1)-1))/N*Fs;
