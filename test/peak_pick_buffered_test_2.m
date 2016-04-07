% File to load in
% Computes peaks of powerspectrum and then performs selective frequency and time
% reassignment.
fname='/tmp/test1.f64';
fi=fopen(fname,'r');
x=fread(fi,Inf,'float64');
x+=randn(length(x),1)*0.01;
fclose(fi);
Fs=16000;
% Window size, for buffering
L=1024;
% Hop size, for buffering
H=256;
% DFT size
N=1024;
X=buffer_ne(x,L,H);
PX=mper(X,'blackman',N,L-1);
newplot(figure(1));
figure(1);
%imagesc((0:(size(PX,2)-1))*H,(0:(N/2-1))/(N)*Fs,a2db(PX(1:(N/2),:)));
%imagesc((0:(size(PX,2)-1)),(0:(N/2-1))/(N)*Fs,a2db(PX(1:(N/2),:)));
hold on
% Cell to store peaks
Pxm=cell();
K=20;
for n=1:K
    Xn=PX(:,n);
    [Xm,Xmi]=lextrem(Xn,'max');
    [Xmin,Xmini]=lextrem(Xn,'min');
    [Xm_,Xmi_]=lmax_parab(p2db(Xn),Xmi);
    [ex_i,Y]=ppvp(Xm_,0.1,'max',0.50,2);
    w=(Xmi_(ex_i)-1)/length(Xn)*(2*pi);
    [X_r,w_r,psi_r,mu_r,t_r,X_r_]=rm_w(X(:,n),w);
    Pxm{n,1}=Xmi_(ex_i);
    Pxm{n,2}=Xm_(ex_i);
    Pxm{n,3}=X_r;
    Pxm{n,4}=w_r;
    Pxm{n,5}=psi_r;
    Pxm{n,6}=mu_r;
    Pxm{n,7}=t_r;
    Pxm{n,8}=X_r_;
    % centre point index calculation
%    xplt=(Xmi_(ex_i)-1);
%    yplt=(Xm_(ex_i))(find(xplt<(N/2)));
    yplt=X_r(find(w_r<pi));
    yplt=a2db(abs(yplt));
    xplt=w_r(find(w_r<pi));
    xplt=xplt/(2*pi)*Fs;
%    xplt=xplt(find(xplt<(N/2)));
%    xplt=xplt/N*Fs;
    xplt_l=length(xplt);
    % Line segment calculation
    t0=-H/2*ones(length(t_r),1);
    t1=H/2*ones(length(t_r),1);
    w0=w_r+psi_r.*t0;
    w1=w_r+psi_r.*t1;
    plot(([t0/H,t1/H]+(n-1)*ones(length(t_r),2))',[w0,w1]'/(2*pi)*Fs);
%    h=scatter((n-1)*ones(xplt_l,1),xplt,5,clip(-1*yplt/50*ones(1,3)));
    h=scatter((n-1)*ones(xplt_l,1),xplt,5,clip(-1*yplt/50*ones(1,3)));
    set(h,'linewidth',2);
end
hold off;
% Example plot
% K=30;
% plot((1:size(PX,1))-1,p2db(PX(:,K)),Pxm{K,1}-1,Pxm{K,2},'o-',Pxm{K,4}/(2*pi)*size(PX,1),Pxm{K,2},'x');
