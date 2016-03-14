% Test the RLS adaptive filter on a signal whose filter varies in time
N1 = 500;
N2 = 500;
N = N1+N2;
x1 = randn(1,N1);
x2 = randn(1,N2);
A = [1];
B1 = [0.5 -0.5];
B2 = [-0.25 0.75];
d1 = filter(B1,A,x1);
d2 = filter(B2,A,x2);
x = [x1 x2];
d = [d1 d2];
p = 2;
lambda=0.97;
[_A,E] = rls(x,d,p,lambda);
t = [1:N] - 1;
plot(t,_A);
