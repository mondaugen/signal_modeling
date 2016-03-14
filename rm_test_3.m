function rm_test_3(path,pmax,N,H,Fs)
% pmax specified in dB Power
newplot(figure(1));
figure(1);
hold on;
[a,phi,w,psi,mu,t,T,F]=read_rm(path,N,H,Fs);
% threshold is percentage of maximum value
thresh=max(max(a))*10^(pmax/20);
for n=1:size(a,2)
    X_a=a(:,n);
    d=diff(X_a);
    % only plot maxima greater than threshold and under nyquist
    idmx=find(([0;(d>0)]&[(d<0);0])&(X_a>thresh)&(w(:,n)<pi));
    X_max=X_a(idmx);
    w_max=w(idmx,n);
    n_max=idmx;
    t_max=t(idmx,n);
    psi_max=psi(idmx,n);
    t0=-H/2*ones(length(idmx),1);
    t1=H/2*ones(length(idmx),1);
    w0=w_max+psi_max.*t0;
    w1=w_max+psi_max.*t1;
    scatter((n-1)*ones(length(idmx),1)*H+t_max,w_max/(2*pi)*Fs,[],idmx);
    plot(([t0+t_max,t1+t_max]+(n-1)*ones(length(idmx),2)*H)',[w0,w1]'/(2*pi)*Fs);
end
hold off;
