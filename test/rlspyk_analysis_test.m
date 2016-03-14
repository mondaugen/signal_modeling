clear;
% Testing a time varying filter
Fs=4410;
f1=1000/Fs;
f2=440/Fs;
r1=0.5;
r2=0.9999;
A1=[1 -2*r1*cos(2*pi*f1) r1*r1];
B1=[(1 - r1*r1)/2 0 -(1 - r1*r1)/2];
A2=[1 -2*r2*cos(2*pi*f2) r2*r2];
B2=[(1 - r2*r2)/2 0 -(1 - r2*r2)/2];
N=44100;
A=interp1([0 N-1],[A1; A2], (0:(N-1)));
B=interp1([0 N-1],[B1; B2], (0:(N-1)));
x=2*rand(1,N)-1;
d=zeros(1,N);
% initial filtering to get state vector
sf=zeros(max(length(A(1,:)),length(B(1,:)))-1,1);
for n=(1:N),
%    [y(n),sf]=filter(B(n,:),A(n,:),x(n),sf);
    [d(n),sf]=filter([1],A(n,:),x(n),sf);
end;
p = 2;
lambda = 0.9999;
[A_,E] = rlspyk_analysis(d,p,lambda);
t = [1:N+1] - 1;
%plot(t,A_,t,[d'; 0]);
plot(t,A_);
