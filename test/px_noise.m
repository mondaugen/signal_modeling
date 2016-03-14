% Compute the power spectra of noise
Fs=44100;
N_sec=20.0;
N=Fs*N_sec;
V1=2;
V2=4;
rho1=sqrt(V1);
rho2=sqrt(V2);
n=(1:N)-1;
t=n/Fs;
T=N/Fs;
x=randn(1,N).*interp1([1 N],[rho1 rho2],(1:N));
M=2048;   % Size of signal chunks
L=512;   % Length of sub-windows
H_w=128;  % Hopsize in Welch power spectrum computation
O_w=L-H_w; % Overlap in Welch power spectrum computation
PX=zeros(M,floor(N/M));
A_est=zeros(1,floor(N/M));
for k=(1:floor(N/M)),
    [PX(:,k),ENBW]=welch(x((k-1)*M+1:k*M),L,O_w,'blackman',M);
end;
figure(1);
imagesc((1:floor(N/M))-1,(0:M-1)/M*Fs,flipud(10*log10(PX)));
P_tot=sum(PX)./ENBW; 
P_noise=P_tot;
figure(2);
plot((1:floor(N/M))-1,[ P_tot(:) P_noise(:) ]);
legend('P_{tot}','P_{noise}');
