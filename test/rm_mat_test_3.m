%fname='~/Documents/sounds/masters_thesis/viola_oboe_violin.raw';
fname='~/Documents/sound/masters_thesis/2016-04-20T18:36:56,777449000-0400_mix.f64le';
fi=fopen(fname,'r');
x=fread(fi,Inf,'float64');
fclose(fi);
L=length(x);
N=2048;
H=512;
P=2;
t=(0:(L-1));
X=buffer_ne(x,N,H);
[Y,w_,psi_,mu_,t_,Y_]=rm(X,'blackman');
Y_plt=cell();
newplot(figure(1));
hold on;
for n=1:size(Y,2)
    S=Y(:,n);
    [lmx,lmxi]=lextrem(abs(S(1:end/2)));
    [ma_i]=lmaxw(lmx,16,1);
    Y_plt{n,1}=lmxi(ma_i);
    Y_plt{n,2}=lmx(ma_i);
    scatter(n*ones(length(ma_i),1),w_(lmxi(ma_i)),[],lmx(ma_i));
end
hold off;

