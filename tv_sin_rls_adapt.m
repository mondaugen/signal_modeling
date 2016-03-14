N=20000;
M=1000;
Fs=1000;
f1=10;
f2=5;
w1=2*pi*f1/Fs;
w2=2*pi*f2/Fs;
amp1=1;
amp2=1;
%x=zeros(N,1);
%x(1)=1;
x=randn(N,1);
%B1=[amp1/2 0 -amp1/2];
B1=[1 0 0];
A1=[1 -2*cos(w1) 1];
%B2=[amp2/2 0 -amp2/2];
B2=[1 0 0];
A2=[1 -2*cos(w2) 1];
B=interp1([1 M],[B1;B2;],1:M);
A=interp1([1 M],[A1;A2;],1:M);
B=[B;interp1([M+1 N],[B2;B2;],(M+1):N);];
A=[A;interp1([M+1 N],[A2;A2;],(M+1):N);];
y=tv_filter(B,A,x);
% carry out linear predicition of next sample from p previous samples using
% linear prediction
p=8;
W=zeros(N,p);
x_=zeros(p,1);
e=zeros(N,1);
c=zeros(N,1);
d=y;
lamb=0.99;
delta=0.001;
w=zeros(p,1);
P=eye(p)/delta;
% predicted value of d
for n=1:N
    c(n)=w'*flipud(x_);
    [w,P,e(n)]=rls_1(d(n),x_,w,P,lamb);
    W(n,:)=w';
    x_=[x_(2:end);y(n)];
end
plot(1:N,y,1:N,c);
