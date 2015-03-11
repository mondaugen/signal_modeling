[d,fs,bps] = ...
 wavread(...
 '~/Documents/sounds/samples/logic-yamaha-grand-piano/047_ped_s_mono.wav');
N = 60000;
d = d(1:N);
x = filter([0 1],[1],d); % delay d by one sample to get input
bet = 0.01;
power = sum(d .^ 2) / N;
p = 50;
mu = bet * 2 / ((p+1)*power);
[_A,E] = lms(x,d,mu,p);
t = [1:N] - 1;
plot(t,_A);
