N=512;
a=1;
phi=0.02;
mu=0.0001;
w=2*pi*10/N;
psi=0.0001;
n_mid=floor(N/2);
t=(1:N)-n_mid-1;
x=a*exp(mu*t+j*(phi+w*t+psi/2*t.^2));
figure(1);
plot(t,real(x));
w=2*pi*(0:(N-1))/N;
[X_,w_,psi_,mu_,t_] = rm(x,'nutall3');
figure(2);
plot((0:(N-1)),[w_/(2*pi) psi_/(2*pi) mu_ abs(X_)]);
legend()
