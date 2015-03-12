N = 10000;
v = randn(1,N);
B = [1];
A = [1 -1.2728 0.81];
d = filter(B,A,v);
p = 2;
lambda = 1.;
x = filter([0 1],[1],d); % delay d by one sample to get input
%[_A,E] = rls(x,d,p,lambda);
[_A,E] = rls_lp_analysis(d,p,lambda);
t = [1:N] - 1;
plot(t,_A);
