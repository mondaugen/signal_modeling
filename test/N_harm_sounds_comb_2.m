% Estimate frequencies by phase unwrapping
% Observe that close frequencies can appear when playing common chords,
% frequencies that are so close they cannot be resolved by simple methods.
Fs=44100;
pchs=[60;64;67];
f0s=midi2hz(pchs);
% amplitudes
a0s=ones(length(pchs),1);
P=floor(Fs/2/min(f0s));
fs=f0s*(1:P);
% partials decay at rate 1/(p^d_rate) where p partial number
d_rate=2;
as=a0s.*1./((1:P).^d_rate);
fs=fs(:)';
as=as(:)';
fs=fs(find(fs<(Fs/2)));
as=as(find(fs<(Fs/2)));
[fs,fsi]=sort(fs);
as=as(fsi);
N=4096;
N_dft=2*N;
H=1024;
a='blackman4';
s_thresh=0.03;
L=N+H;
n=(0:(N-1))';
n_dft=(0:(N_dft-1));
l=(0:(L-1))';
x=exp(j*2*pi*(l/Fs*fs)).*as;
x=sum(x,2);
% low-pass filter for noise
fc=5000;
K_lp=fc*2*pi/Fs;
a_flt=K_lp/(K_lp+1);
P_no=-Inf;
x_no=randn(L,1)*sqrt(10^(P_no/20)*length(fs));
y_no=filter([a_flt],[1,-(1-a_flt)],x_no);
x+=y_no;
%w=sum_cos_win_t(n,[0.42,0.5,0.08],4096);
w=sum_cos_win_t(n,a,N);
x_0=x(1:N);
x_1=x(H:(N+H-1));
X_0=fft(x_0.*w,N_dft)/sum(w);
X_1=fft(x_1.*w,N_dft)/sum(w);
[Xm,Xmi]=lmax(abs(X_0));
[Xm_,Xmi_,p_]=lmax_parab(log(abs(X_0)),Xmi);
Xmi_=Xmi_-1;
figure(1);
plot(n_dft,log(abs(X_0)),Xmi-1,log(Xm),'o',Xmi_,Xm_,'o');
plot_vert_lines(fs,min(log(abs(X_0))),max(log(abs(X_0))));
% using estimated frequencies check correlation with frequency domain window
W=zeros(N_dft,length(Xmi_));
W0=sum_cos_win_f(0,a,N);
% compute the sidelobe amplitude (for blackman, around the 10th bin is where the
% sidelobe starts)
Wthresh=sum_cos_win_f(10/N_dft*2*pi,a,N)/W0;
for k=1:length(Xmi_)
    W(:,k)=sum_cos_win_f((n_dft(:)/N_dft-Xmi_(k)/N_dft)*2*pi,a,N)/W0;
end
W=W.*(abs(W)>Wthresh);
W=W.';
s=(abs(W)*abs(X_0))./(norm(abs(W),2,"rows")*norm(abs(X_0),2));
%s=(abs(W)*abs(X_0))
%s=W*X_0./(norm(W,2,"rows")*norm(X_0,2));
%s=W*X_0;
Xmi_threshed=Xmi_(find(s>s_thresh));
Xm_threshed=Xm_(find(s>s_thresh));
figure(2);
plot(n_dft/N_dft*Fs,log(abs(X_0)),Xmi_threshed/N_dft*Fs,Xm_threshed,'o');
plot_vert_lines(fs,0,max(s));
figure(3);
stem(Xmi_/N_dft*Fs,log(s),"basevalue",min(log(s)));
%[X_0_r,w_0_r]=rm(x_0,'hanning');
%plot(n/N*44100,abs(X_0)./sum(w),w_0_r/(2*pi)*44100,abs(X_0_r),'o',fs,ones(length(fs),1),'x');
%plot(n,abs(X_0)./sum(w))
