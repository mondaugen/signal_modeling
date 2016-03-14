clear;
% Testing a time varying filter
Fs=44100;
f1=1000/Fs;
f2=440/Fs;
r1=0.5;
r2=0.9999;
A1=[1 -2*r1*cos(2*pi*f1) r1*r1];
B1=[(1 - r1*r1)/2 0 -(1 - r1*r1)/2];
A2=[1 -2*r2*cos(2*pi*f2) r2*r2];
B2=[(1 - r2*r2)/2 0 -(1 - r2*r2)/2];
N=4410;
A=interp1([0 N-1],[A1; A2], (0:(N-1)));
B=interp1([0 N-1],[B1; B2], (0:(N-1)));
x=2*rand(1,N)-1;
d=zeros(1,N);
% initial filtering to get state vector
sf=zeros(max(length(A(1,:)),length(B(1,:)))-1,1);
for n=(1:N),
    [d(n),sf]=filter(B(n,:),A(n,:),x(n),sf);
%    [d(n),sf]=filter([1],A(n,:),x(n),sf);
end;
lambda = 0.9999;
p = 2;
% L, the index of the sample to insert
L = 3000;
% K, the number of times to insert new samples
K = 1400;
[A_,E] = rlspyk_analysis(d,p,lambda);
t = [1:N+1] - 1;
% Resynthesize original residual
[B_,D] = rlspyk_synthesis(E,p,lambda);
for k=(1:K),
    % Calculate sample to insert into residual
    e_=D(L)-B_(L+1,1)*B_(L,:)*D((L-1)-(0:(p-1)))-B_(L+1,(2:end))*D(L-(1:(p-1)))-E(L);
    e_=e_/B_(L+1,1);
    % Insert sample
    E=[E(1:(L-1)) e_ E((L:end))];
    % increment L
    L=L+1;
end;
% Try synthesizing again
[B__,d] = rlspyk_synthesis(E,p,lambda);
t_ = (1:size(B__,1))-1;
figure(2);
plot(t_(2:end),B__(2:end,:),t_(2:end),d);
