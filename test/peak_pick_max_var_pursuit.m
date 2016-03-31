% Generate random signal
clear;
N=1000;
Y=randn(N,1);
X=Y;
ex_states=mark_extrem(X);
figure(1);
plot(1:N,X,1:N,ex_states,'o');
ex_i=(1:N);
va=abs(X-[X(1);X(1:(N-1))])+abs(X-[X(2:N);X(N)]);
tot_va=sum(va);
thresh=tot_va*0.01;
t=1;
fig=figure(2);
set(fig, 'visible', 'off');
while (tot_va > thresh)
    mi_i=find(ex_states==-1);
    [va_ma,va_ma_i]=max(va(mi_i));
    va_ma_i=ex_i(mi_i(va_ma_i));
    ex_i=ex_i(find(ex_i != va_ma_i));
    X=Y(ex_i);
    ex_states=mark_extrem(X);
    N-=1;
    va=abs(X-[X(1);X(1:(N-1))])+abs(X-[X(2:N);X(N)]);
    tot_va=sum(va);
%    pause;
%    sleep(1);
    filename=sprintf('/tmp/vp_%d.png',t);
    t+=1;
    plot(1:length(Y),Y,ex_i,X,'o-');
%    print(filename);
    saveas(2,filename,'png');
end
