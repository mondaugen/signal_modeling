f=400;
fs=16000;
T=1;
t=0:(1/fs):1;
H=16;
x=cos(2*pi*f*t);
[S,T,F]=stft_ne(x,1024,H,'hanning',512,16000);
plot(0:H:(length(x)-1),arg(S(14,:)))
