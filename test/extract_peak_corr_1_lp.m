% Use the correlation of the FT of two signals to extract a sine from noise
clear;
Fs=16000;
f=3000;
N=1024;
N_dft=8192;
n=(0:(N-1));
n_dft=(0:(N_dft-1))/N_dft*N;
n_dft=n_dft(:);
n=n(:);
P_sin=-30;
a_sin=10^(P_sin/20);
P_no=0;
x=zeros(N,1);
x_no=randn(N,1)*sqrt(10^(P_no/20));
% low-pass filter for noise
fc=500;
K_lp=fc*2*pi/Fs;
a_flt=K_lp/(K_lp+1);
y_no=filter([a_flt],[1,-(1-a_flt)],x_no);
x+=y_no;
x+=a_sin*exp(j*2*pi*f*n/Fs);
%x+=cos(2*pi*f*n/Fs);
w=sum_cos_win_t(n,'blackman4-min',N);
[X,ENBW]=mper(x,w,N_dft,N);
[Xm,Xmi]=lmax(X);
[Xm_,Xmi_]=lmax_parab(p2db(X),Xmi);
[W_]=mper(exp(j*2*pi*f*n/Fs),w,N_dft,N);
B=4;
Wc=zeros(N_dft,length(Xmi_));
for k_=1:length(Xmi_)
    l_=(Xmi_(k_)-1);
    Wc(:,k_)=db2p(Xm_(k_))*mper(exp(j*2*pi*l_/N_dft*n),w,N_dft,N);
    Wc(:,k_).*=((n_dft >= (l_/N_dft*N-B)) & (n_dft <= (l_/N_dft*N+B)));
end
figure(1);
plot(n_dft/N*Fs,p2db(X))
plot_vert_lines(f,min(p2db(X)),max(p2db(X)));
Wc=Wc.';
s=Wc*X./(norm(Wc,2,'rows')*norm(X,2));
figure(2);
plot((Xmi_-1)/N_dft*Fs,s);
plot_vert_lines(f,min(s),max(s));
