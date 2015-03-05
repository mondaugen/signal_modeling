% compares various prony model orders of an arbitrary signal
t0 = 100;
N = 200;
t = (t0:(t0+N-1));
Pmin = 20; % Min number of poles
Pmax = 100; % Max number of poles
p = [Pmin:Pmax]; % number of poles
q = 1; % number of zeros
fpath = ...
'~/Documents/sounds/samples/logic-yamaha-grand-piano/047_ped_s_mono-8k.wav';
[x,fs,bps] = wavread(fpath);
x = x(t+1);
plot(t,x,'b');
hold on;
for _p=p,
    [y,err] = prony_impulse_model(x,_p,q);
    plot(t,y,'g');
end;
hold off;
