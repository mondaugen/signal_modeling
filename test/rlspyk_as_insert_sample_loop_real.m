clear;
[d,fs]=wavread('tmp/maybe_if.wav');
% left channel only
d=d(:,1);
N=length(d);
lambda = 0.9999;
p = 50;
% L, the index of the sample to insert
L = 200;
% K, the number of times to insert new samples
K = 10;
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
wavwrite(d,fs,'tmp/maybe_if_extended.wav');
