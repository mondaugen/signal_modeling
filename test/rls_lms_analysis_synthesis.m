N = 10000;
v = randn(1,N);
B = [1];
A = [1 -1.2728 0.81];
d = filter(B,A,v);
lambda = 0.9999;
p_rls = 4;
p_lms = 60;
be = 0.001;
t = [1:N] - 1;

[A_,U] = rlspyk_analysis(d,p_rls,lambda);
subplot(2,2,1);
plot(t,A_(2:end,:));

'power reduction of signal'
6*log(sum(U .^ 2)/sum(d .^ 2))/log(2)

[A_,E] = nlms_lp_analysis(U,be,p_lms);
subplot(2,2,3);
plot(t,A_);

'power reduction of signal'
6*log(sum(E .^ 2)/sum(d .^ 2))/log(2)

[A_,U] = nlms_lp_synthesis(E,be,p_lms);
subplot(2,2,4);
plot(t,A_);

[A_,D] = rlspyk_synthesis(U,p_rls,lambda);
subplot(2,2,2);
plot(t,A_(2:end,:));