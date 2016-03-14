N_ARGS = 5;
args=argv();
if (nargin != N_ARGS)
    fprintf(stderr,...
    ["usage: " program_name() " N H W L format\n" ...
     "N : FFT size\n" ...
     "H : hop size\n" ...
     "W : window type\n" ...
     "    can be 'hanning', 'hamming'," ...
     " 'blackman', 'bartlett', or 'rectangular'\n" ...
     "L : window length\n" ...
     "format : the format of the output\n" ...
     "         can be\n" ...
     "         'raw' (doubles) or\n" ...
     "         'human' (printf doubles of magnitude spectrum)\n\n" ...
     "Reads stream of values of type \"float64\" from stdin."]);
    exit(-1);
end
N=str2num(args{1});
H=str2num(args{2});
W=args{3};
L=str2num(args{4});
x=zeros(L,1);
format=args{5};
if (L > N)
    fprintf(stderr,'Window sizes greater than N not yet supported.');
    exit(-1);
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
end
w=w(1:L);
W_0=sum(w);
% Space to interleave real and imaginary parts
X_=zeros(2*N,1);
% Load in enough samples to fill half of x so that the first frame's centre is
% aligned with the first sample
[tmp,count]=fread(stdin,L/2,'double');
tmp=tmp(:);
x((L/2+1):(L/2+count))=tmp;
if (count != L/2)
    x((L/2+1+count):L)=zeros(L/2-count,1);
end
do
        X=fft(x.*w,N)/W_0;
    switch (format)
        case ('human')
            printf("%f ",abs(X));
            printf("\n");
        case ('raw')
            X_(1:2:end)=real(X);
            X_(2:2:end)=imag(X);
            fwrite(stdout,X_,'double');
    end
        % shift over values
        x(1:(L-H))=x((H+1):L);
        [tmp,count]=fread(stdin,H,'double');
        tmp=tmp(:);
        x((L-H+1):(L-H+count))=tmp;
until (count != H)
% pad with zeros until buffer full
x((L-H+1+count):L)=zeros(H-count,1);
X=fft(x.*w,N)/W_0;
switch (format)
    case ('human')
        printf("%f ",abs(X));
        printf("\n");
    case ('raw')
        X_(1:2:end)=real(X);
        X_(2:2:end)=imag(X);
        fwrite(stdout,X_,'double');
end

