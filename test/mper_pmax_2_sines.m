N=512;
N_fft=2048;
A1=0.5;
rho=sqrt((1/2)*A1^2); % SNR is equal to unity
f1=100.1;
f2=200.2;
n=(0:(N-1));
K=100;
Px=zeros(1,N_fft);
for k=(1:K),
    phi1=rand(1,1)*2*pi-pi;
    phi2=rand(1,1)*2*pi-pi;
    x=A1*cos((n/N*2*pi*f1.+phi1));
    %x+=cos((n/N*2*pi*f2.+phi2));
    x+=randn(1,N)*rho;
    [Px_,ENBW]=mper(x,'blackman',N_fft);
    Px=Px.+Px_';
end;
% take average
Px=Px/K;
Px_log=10*log10(Px);
[x_max,y_max]=parab_max(Px_log);
n_fft=(1:N_fft)-1;
plot(n_fft,Px_log,x_max,y_max);
disp('total power')
sum(Px)/ENBW
