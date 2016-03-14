L=256;
N=512;
N_fft=N;%2048;
O=128;
A1=0.5;
rho=sqrt((1/2)*A1^2); % SNR is equal to unity
f1=100.1;
f2=200.2;
n=(0:(N-1));
K=100;
Px=zeros(1,N_fft);
phi1=rand(1,1)*2*pi-pi;
%phi2=rand(1,1)*2*pi-pi;
x=A1*cos((n/N*2*pi*f1.+phi1));
%x+=cos((n/N*2*pi*f2.+phi2));
x+=randn(1,N)*rho;
[Px,ENBW]=welch(x,L,O,'blackman');
Px_log=10*log10(Px);
[x_max,y_max]=parab_max(Px_log);
n_fft=(1:N_fft)-1;
plot(n_fft,Px_log,x_max,y_max);
disp('total power')
P_tot=sum(Px)/ENBW
disp('estimated magnitude of sinusoid')
A_est=sqrt(4*10^(y_max/10))
disp('estimated power of sinusoid')
P_A=0.5*A_est^2
disp('estimated power of noise')
P_noise=P_tot-P_A
disp('harmonic SNR')
P_A/P_noise
