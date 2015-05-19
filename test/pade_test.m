% model 1000 samples with bi-quad using pade's approx
N = 1000;
t = (1:N) - 1;
p = 2; % number of poles
q = 2; % number of zeros
fpath = ...
'~/Documents/sounds/samples/logic-yamaha-grand-piano/047_ped_s_mono.wav';
[x,fs,bps] = wavread(fpath);
x = x(1:N);
[a,b] = pade(x,p,q);
d = zeros(1,N);
d(1)=1;
y = filter(b,a,d);
plot(t,x,y);
