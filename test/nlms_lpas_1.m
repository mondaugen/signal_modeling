[d,fs,bps] = ...
 wavread(...
 '~/Documents/sounds/samples/logic-yamaha-grand-piano/047_ped_s_mono.wav');
N = 20000;
d = d(1:N);
bet = 0.001;
p = 50;
[_A,E] = nlms_lp_analysis(d,bet,p);
t = [1:N] - 1;
subplot(2,1,1);
plot(t,_A);
[_B,D] = nlms_lp_synthesis(E,bet,p);
subplot(2,1,2);
plot(t,_B);
