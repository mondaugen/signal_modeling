clear;
% Sample rate
Fs=8000;
N = 10000;
v = randn(1,N);
B = [1];
A = [1 -1.2728 0.81];
d = filter(B,A,v);
lambda = 0.9999;
p = 2;
% L, the index of the sample to insert
L = 5000;
% K, the number of times to insert new samples
K = 2000;
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
    L=L-2;
end;
% Try synthesizing again
[B__,d] = rlspyk_synthesis(E,p,lambda);
t_ = (1:size(B__,1))-1;
figure(2);
% Normalize 
d=d./max([max(d) -min(d)]);
plot(t_(2:end),B__(2:end,:),t_(2:end),d);
wavwrite(d,Fs,'tmp/filled_noise.wav');
