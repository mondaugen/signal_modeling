% model 1000 samples with bi-quad using pade's approx
t0 = 1000;
N = 1000;
t = (t0:(t0+N-1));
p = 20; % number of poles
q = 1; % number of zeros
fpath = ...
'~/Documents/sounds/samples/logic-yamaha-grand-piano/047_ped_s_mono.wav';
[x,fs,bps] = wavread(fpath);
x = x(t+1);
[a,b,err] = prony(x,p,q);
d = zeros(1,N);
d(1)=1;
y = filter(b,a,d);
plot(t,x,'b',t,y,'g');
