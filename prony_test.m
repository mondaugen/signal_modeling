function prony_test(fpath,p,q,t0,N);
% model 1000 samples with bi-quad using prony's method
% p is number of poles
% q is number of zeros
if nargin<5, N=1000; end;
if nargin<4, t0=1000; end;
t = (t0:(t0+N-1));
[x,fs,bps] = wavread(fpath);
if size(x)(2) > 1,
    printf('Multichannel file given. Only reading 1st channel.\n');
end;
x = x(:,1);
x = x(t+1);
[a,b,err] = prony(x,p,q);
d = zeros(1,N);
d(1)=1;
y = filter(b,a,d);
plot(t,x,'b',t,y,'g');
