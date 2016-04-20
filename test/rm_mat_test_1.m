N=2048;
H=512;
P=2;
L=100000;
t=(0:(L-1));
a=unifrnd(0.1,1,[P,1]);
phi=unifrnd(-pi,pi,[P 1]);
w=unifrnd(0,pi,[P 1]);
mu=unifrnd(log(1e-1./a)/L,log(1./a)/L);
psi=unifrnd(-w./L,(pi-w)/L);
x=sum(a.*exp(mu*t+j*(phi+w*t+psi/2*t.^2)));
X=buffer_ne(x,N,H);
[Y,w_,psi_,mu_,t_,Y_]=rm(X);
w_plt=w_(:);
t_plt=ones(N,1)*((1:size(Y,2))-1)*H+t_;
t_plt=t_plt(:);
i_plt=find(abs(Y_)>0.0001);
h=scatter(t_plt(i_plt),w_plt(i_plt),[],'k');
set(h,'linewidth',1);