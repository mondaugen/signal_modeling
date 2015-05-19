% Test if periodogram needs to be scaled by N or not.
N=256;
f=100.5;
n=(0:(N-1));
K=100;
Px=zeros(1,N);
for k=(1:K),
    phi=rand(1,1)*2*pi-pi;
    x=sin(n/N*2*pi*f.+phi);%+randn(1,N);
    % Divide again by N because of integral of the window (?)
    Px=Px.+mper(x,'boxcar')';
end;
Px=Px/K;
a=sqrt(2*interp1((1:(N)),Px,f+1,'pchip'))
b=sqrt(2*interp1((1:(N)),Px,N-f+1,'pchip'))
sqrt(sum(a^2+b^2))
stem(n,Px);
