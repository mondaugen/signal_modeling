function [x,S,S2,n_max,X_max,w_max,t_max]=rm_test_2(thresh=0.5);
    N=2048;
    H=512;
    P=2;
    L=10000;
    t=(0:(L-1));
    a=unifrnd(0.1,1,[P,1]);
    phi=unifrnd(-pi,pi,[P 1]);
    w=unifrnd(0,pi,[P 1]);
    mu=unifrnd(log(1e-1./a)/L,log(1./a)/L);
    %psi=unifrnd(-0.5*pi/L,0.5*pi/L,[P 1]);
    psi=unifrnd(-w./L,(pi-w)/L);
    x=sum(a.*exp(mu*t+j*(phi+w*t+psi/2*t.^2)));
    [S,T]=stft_ne(x,N,H,'hanning',N);
    w_plt=2*pi*(0:(N-1))/N;
    w_plt=w_plt(:);
    newplot(figure(1));
    figure(1);
    imagesc(T,w_plt,abs(S));
    S2=fproc_eo(x,N,H,'rectangular',N,@_f,N);
    X_=S2{1,1};
    w_=S2{2,1};
    t_=S2{5,1};
    psi_=S2{3,1};
    n_max=cell(1,size(X_,2));
    X_max=cell(1,size(X_,2));
    w_max=cell(1,size(X_,2));
    t_max=cell(1,size(X_,2));
    psi_max=cell(1,size(X_,2));
    newplot(figure(2));
    figure(2);
    hold on;
    for n=1:size(X_,2)
        X=X_(:,n);
        w__=w_(:,n);
        t__=t_(:,n);
        psi__=psi_(:,n);
        X_a=abs(X);
        d=diff(abs(X));
        idmx=find(([0;(d>0)]&[(d<0);0])&(X_a>thresh));
        X_max(1,n)=abs(X(idmx));
        w_max(1,n)=w__(idmx);
        n_max(1,n)=idmx;
        t_max(1,n)=t__(idmx);
        psi_max(1,n)=psi__(idmx);
        psi_max__=psi__(idmx);
        w_max__=w__(idmx);
        t_max__=t__(idmx);
        t__0=-H/2*ones(length(idmx),1);
        t__1=H/2*ones(length(idmx),1);
%        w__0=max(min(w_max__+psi_max__.*t__0,w_plt(idmx)+1/N*2*pi),w_plt(idmx)-1/N*2*pi);
%        w__1=max(min(w_max__+psi_max__.*t__1,w_plt(idmx)+1/N*2*pi),w_plt(idmx)-1/N*2*pi);
        w__0=w_max__+psi_max__.*t__0;
        w__1=w_max__+psi_max__.*t__1;
        scatter((n-1)*ones(length(idmx),1)*H+t_max{1,n},w_max{1,n},[],idmx);%X_max{1,n});
        plot(([t__0+t_max{1,n},t__1+t_max{1,n}]+(n-1)*ones(length(idmx),2)*H)',[w__0,w__1]');
%        scatter((n-1)*ones(length(idmx),1)*H,w_max{1,n},[],'r');%X_max{1,n});
    end
    hold off;
end

function [X_,w_,psi_,mu_,t_] = _f(x,n)
    [X_,w_,psi_,mu_,t_]=rm(x,'hanning');
end
