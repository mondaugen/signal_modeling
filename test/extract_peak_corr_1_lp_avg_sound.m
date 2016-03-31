% Extract sines from real sound
clear;
N=4096;
%filepath='~/Documents/sound/masters_thesis/FhrnA#2_2.27sec.wav';
filepath='~/Documents/sound/masters_thesis/VlnC6_3.88sec.wav';
%filepath='~/Documents/sound/masters_thesis/ELECTRIC_BASS_DEEP_E3.wav';
[x,Fs]=wavread(filepath);
%f=116.5
f=1046.5;
%f=82.407
fs=f*(1:floor((Fs/2)/f));
n0=Fs;
x=x(n0:(N+n0-1),1);
f=3000;
N_dft=8192;
n=(0:(N-1));
n_dft=(0:(N_dft-1))/N_dft*N;
n_dft=n_dft(:);
n=n(:);
P_sin=-36;
a_sin=10^(P_sin/20);
P_no=0;
%x+=cos(2*pi*f*n/Fs);
w=sum_cos_win_t(n,'blackman',N);
%X=zeros(N_dft,1);
[X,ENBW]=mper(x,w,N_dft,N);
[Xm,Xmi]=lextrem(X,'max');
[Xmin,Xmini]=lextrem(X,'min');
[Xm_,Xmi_]=lmax_parab(p2db(X),Xmi);
[W_]=mper(exp(j*2*pi*f*n/Fs),w,N_dft,N);
B=3;
Wc=zeros(N_dft,length(Xmi_));
for k_=1:length(Xmi_)
    l_=(Xmi_(k_)-1);
    Wc(:,k_)=mper(exp(j*2*pi*l_/N_dft*n),w,N_dft,N); %db2p(Xm_(k_))*
    Wc(:,k_).*=((n_dft >= (l_/N_dft*N-B)) & (n_dft <= (l_/N_dft*N+B)));
end
%for k_=1:length(n_dft)
%    l_=k_-1;
%    Wc(:,k_)=mper(exp(j*2*pi*l_/N_dft*n),w,N_dft,N); %db2p(Xm_(k_))*
%    Wc(:,k_).*=((n_dft >= (l_/N_dft*N-B)) & (n_dft <= (l_/N_dft*N+B)));
%end
figure(1);
plot(n_dft/N*Fs,p2db(X),(Xmi_-1)/N_dft*Fs,Xm_,'o-')
plot_vert_lines(fs,min(p2db(X)),max(p2db(X)));
Wc=Wc.';
s=Wc*X./(norm(Wc,2,'rows')*norm(X,2));
figure(2);
plot((Xmi_-1)/N_dft*Fs,p2db(s));
plot_vert_lines(fs,min(p2db(s)),max(p2db(s)));
figure(3);
x_smooth=inv_var_filt(p2db(X),2,1);
plot(n_dft/N*Fs,x_smooth);
