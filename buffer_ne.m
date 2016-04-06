function [S,T]=buffer_ne(s,N,H)
% [S]=buffer_ne(x,N,H)

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
Fs=1;
if (L_end < N/2)
    x((N/2+1):(N/2+L_end))=s;
    x((N/2+1+L_end):N)=zeros(N/2-L_end,1);
else
    x((N/2+1):N)=s(1:N/2);
end
n=n+H;
while (n <= length(s))
    T=[T (n-1)/Fs];
    S=[S x];
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
S=[S x];
