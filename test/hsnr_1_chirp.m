% Compute the harmonic signal-to-noise-ratio
Fs=44100;
f1_0=100;
f1_1=10000;
A1=0.5;
phi0=rand(1,1)*2*pi-pi;
N_sec=2.0;
N=Fs*N_sec;
SNR=0;
rho=sqrt((A1^2)/(2*10^(SNR/10)));
n=(1:N)-1;
t=n/Fs;
T=N/Fs;
phi_n=phi0+2*pi*(f1_0*t+(f1_1-f1_0)*t.^2/(2*T));
%phi_n=phi0+2*pi*(f1_0+(f1_1-f1_0)*t/T);
x=A1*cos(phi_n)+randn(1,N)*rho;
M=2048;   % Size of signal chunks
L=512;   % Length of sub-windows
H_w=128;  % Hopsize in Welch power spectrum computation
O_w=L-H_w; % Overlap in Welch power spectrum computation
PX=zeros(M,floor(N/M));
A_est=zeros(1,floor(N/M));
for k=(1:floor(N/M)),
    [PX(:,k),ENBW]=welch(x((k-1)*M+1:k*M),L,O_w,'blackman',M);
    [x_m,y_m]=parab_max(10*log10(PX(:,k)));
    A_est(k)=sqrt(4*10^(y_m/10));
end;
figure(1);
imagesc((1:floor(N/M))-1,(0:M-1)/M*Fs,flipud(10*log10(PX)));
P_tot=sum(PX)./ENBW; 
P_A=0.5*A_est.^2;
P_noise=P_tot-P_A;
HSNR=P_A./P_noise;
figure(2);
plot((1:floor(N/M))-1,[P_tot(:) A_est(:) P_A(:) P_noise(:) 10*log10(HSNR(:))]);
legend('P_{tot}','A_{est}','P_{A}','P_{noise}','HSNR');
