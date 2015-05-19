clear;
N = 10000;
v = randn(1,N);
B = [1];
A = [1 -1.2728 0.81];
d = filter(B,A,v);
lambda = 0.9999;
p = 2;
% L, the index of the sample to insert
L = 1000;
[A_,E] = rlspyk_analysis(d,p,lambda);
t = [1:N+1] - 1;
figure(1);
plot(t,A_);
% Resynthesize original residual
[B_,D] = rlspyk_synthesis(E,p,lambda);
% Calculate sample to insert into residual
e_=D(L)-B_(L+1,1)*B_(L,:)*D((L-1)-(0:(p-1)))-B_(L+1,(2:end))*D(L-(1:(p-1)))-E(L);
e_=e_/B_(L+1,1);
size(E)
% Insert sample
E=[E(1:(L-1)) e_ E((L:end))];
% Try synthesizing again
[B__,D_] = rlspyk_synthesis(E,p,lambda);
t_ = [1:N+2] - 1;
figure(2);
plot(t_,B__);

