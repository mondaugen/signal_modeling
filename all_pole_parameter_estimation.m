Fs=44100;
f=500;
T60=2.0;
r_T60=10^(-3);
N=Fs*T60;
r=r_T60^(1/N);
phi=2*pi*f/Fs;
w=r*exp(j*phi);
%b=[1 - 0.5*(w+w')];
%a=[1 -(w+w') (w*w')];
%G=abs(sum(b.*((w/r).^[0 -1])))/abs(sum(a.*((w/r).^[0 -1 -2])));
b=[1 -r*cos(phi)];
a=[1 -2*r*cos(phi) r^2];
figure(1)
freqz(b,a,8192);
x=[1 zeros(1,N-1)];
y=filter(b,a,x);
figure(2)
plot(1:length(y),y);
[a_,b_]=prony(y(1:100),2,1);
y_=filter(b_,a_,x);
figure(3)
plot(1:length(y_),y_);
