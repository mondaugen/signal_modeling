% File to load in
fname='/tmp/c.f64';
fi=fopen(fname,'r');
x=fread(fi,Inf,'float64');
x+=randn(length(x),1)*0.1;
fclose(fi);
Fs=16000;
% Window size, for buffering
L=1024;
% Hop size, for buffering
H=256;
X=buffer_ne(x,L,H);
PX=mper(X,'blackman',2048,L-1);
newplot(figure(1));
figure(1);
imagesc(a2db(PX));
hold on
% Cell to store peaks
Pxm=cell();
for n=1:size(PX,2)
    Xn=PX(:,n);
    [Xm,Xmi]=lextrem(Xn,'max');
    [Xmin,Xmini]=lextrem(Xn,'min');
    [Xm_,Xmi_]=lmax_parab(p2db(Xn),Xmi);
    [ex_i,Y]=ppvp(Xm_,0.1,'max',0.50,2);
    Pxm{n,1}=Xmi_(ex_i);
    Pxm{n,2}=Xm_(ex_i);
    scatter(n*ones(length(ex_i),1),Xmi_(ex_i),5);
end
hold off;
