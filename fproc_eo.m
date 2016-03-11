function [S]=fproc_eo(s,N,H,W,L,f,Q)
% Frame by frame processing.
%
% f is function that accepts frame of s and time at start of frame to be
% processed. f can return any number of outputs, but the size of each output
% should not change between calls.
% Q is the size of each output and is used to preallocate memory for speed.
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
        w=ones(L+1,1);
    otherwise
        error(sprintf('Bad window: %s',W));
end
w=w(1:L);
% pad with zeros
w=[w;zeros(N-L,1)];
W_0=sum(w);
% space to load in samples
x=zeros(N,1);
% space to write output
M=nargout(f);
L_tot=length(s);
L_frames=L_tot/(N-H)+1;
S=cell(M,1);
for m=1:M
    S{m,1}=zeros(Q,L_frames);
end
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
n_frame=1;
while (n <= length(s))
    X=nthargout([1:M],f,x,n);
    for m=1:M
        S{m,1}(:,n_frame)=X{1,m};
    end
    % shift over values
    x(1:(N-H))=x((H+1):N);
    L_end=length(s)-n+1;
    if (L_end < H)
        x((N-H+1):(N-H+L_end))=s(n:end);
    else
        x((N-H+1):N)=s(n:(n+H-1));
    end
    n=n+H;
    n_frame=n_frame+1;
end
% pad with zeros until buffer full
x((N-H+1+L_end):N)=zeros(H-L_end,1);
X=nthargout([1:M],f,x,n);
for m=1:M
    S{m,1}(:,n_frame)=X{1,m};
end
