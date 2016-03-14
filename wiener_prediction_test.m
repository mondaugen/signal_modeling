clear;
w=0.1*2*pi;
a=0.99;
A=[1 -2*a*cos(w) a^2];
B=[1];
N=1000;
x=randn(N,1);
y=filter(B,A,x);
p=2;
R=acorrm(y,512)(1:p,1:p);
r=acorr(y,512)(2:(p+1));
w=linsolve(R,r);
x_=filter([1 -w'],B,y);
Rx_=acorr(x,128);
