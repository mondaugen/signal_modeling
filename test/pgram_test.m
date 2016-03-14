% Test if periodogram needs to be scaled by N or not.
N=4096;
f=100;
n=(0:(N-1));
K=1;
Px=zeros(1,N);
for k=(1:K),
    phi=rand(1,1)*2*pi-pi;
    x=sin(n/N*2*pi*f.+phi)+randn(1,N);
    % Divide again by N because of integral of the window (?)
    Px=Px.+(abs(fft(x)).^2)/(N.^2);
end;
Px=Px/K;
Px(f+1)
Px(N-f+1)
