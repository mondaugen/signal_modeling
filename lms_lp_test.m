N = 10000;
v = randn(1,N);
B = [1];
A = [1 -1.2728 0.81];
d = filter(B,A,v);
x = filter([0 1],[1],d); % delay d by one sample to get input
bet = 0.01;
power = sum(d .^ 2) / N;
p = 2;
mu = bet * 2 / ((p+1)*power);
[_A,E] = lms(x,d,mu,p);
t = [1:N] - 1;
plot(t,_A(:,1),_A(:,2));
