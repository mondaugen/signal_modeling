% Estimate frequencies by phase unwrapping
% Observe that close frequencies can appear when playing common chords,
% frequencies that are so close they cannot be resolved by simple methods.
Fs=44100;
pchs=[60;64;67];
f0s=midi2hz(pchs);
P=floor(Fs/2/min(f0s));
fs=f0s*(1:P);
fs=fs(:)';
fs=fs(find(fs<(Fs/2)));
fs=sort(fs);
N=4096;
H=1024;
L=N+H;
n=(0:(N-1))';
l=(0:(L-1))';
x=exp(j*2*pi*(l/Fs*fs));
x=sum(x,2);
w=sum_cos_win_t(n,[0.42,0.5,0.08],4096);
x_0=x(1:N);
x_1=x(H:(N+H-1));
X_0=fft(x_0.*w);
X_1=fft(x_1.*w);
Xd=X_1./X_0;
F_Xd=arg(Xd);
for k=2:(length(F_Xd)/2)
    while(F_Xd(k) < (2*pi*H/N*(k-1.5)))
        F_Xd(k)+=2*pi;
    end
end
F_Xd=F_Xd/(H*2*pi)*Fs;
[sF_Xd,siF_Xd]=sort(F_Xd);
plot(F_Xd(siF_Xd),abs(X_0)(siF_Xd)/sum(w));
plot_vert_lines(fs,0,2);
%[X_0_r,w_0_r]=rm(x_0,'hanning');
%plot(n/N*44100,abs(X_0)./sum(w),w_0_r/(2*pi)*44100,abs(X_0_r),'o',fs,ones(length(fs),1),'x');
%plot(n,abs(X_0)./sum(w))
