% See Hayes p. 152
N=80;
n_d=5;
n=(0:(N-1));
%h=sin((n - n_d)*pi/2)./((n-n_d)*pi);
h=sinc((n-n_d)/2)*2;
p=5;
q=5;
[a,b,err]=prony(h(1:11),p,q);
freqz(b,a);
