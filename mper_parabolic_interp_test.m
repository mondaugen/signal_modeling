N=512;
N_fft=2048;
f=100.125;
n=(0:(N-1));
K=8;
Px=zeros(1,N_fft);
for k=(1:K),
    phi=rand(1,1)*2*pi-pi;
    x=cos((n/N*2*pi*f.+phi))+randn(1,N);
    [Px_,ENBW]=mper(x,'blackman',N_fft);
    Px=Px.+Px_';
end;
% take average
Px=Px/K;
Px_log = 10*log10(Px);
[ma,ind]=max(Px_log)
% interpolate to estimate true maximum
indices=mod((ind-1:ind+1)-1,N_fft)+1;
p_ma_est_log_poly=polyfit(indices,Px_log(indices),2)
p_ma_est_log_poly_dx=polyder(p_ma_est_log_poly);
p_ma_est_log_x=-p_ma_est_log_poly_dx(2)/(p_ma_est_log_poly_dx(1))
p_ma_est_log_y=polyval(p_ma_est_log_poly,p_ma_est_log_x)
p_ma_est=10^(p_ma_est_log_y/10)
%a=sqrt(2*interp1((1:(N)),Px,f+1,'pchip'))
%b=sqrt(2*interp1((1:(N)),Px,N-f+1,'pchip'))
%sqrt(sum(a^2+b^2))
n_fft=(0:(N_fft-1));
plot(n_fft+1,Px_log,(ind-1:0.01:ind+1),polyval(p_ma_est_log_poly,(ind-1:0.01:ind+1)));
disp('total power')
sum(Px)/ENBW

