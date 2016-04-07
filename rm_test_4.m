function rm_test_4(path,pmax,N,H,Fs)
% pmax specified in dB Power
newplot(figure(1));
figure(1);
hold on;
[a,phi,w,psi,mu,t,T,F]=read_rm([path,'.rm'],N,H,Fs);
% threshold is percentage of maximum value
thresh=max(max(a))*10^(pmax/20);
k_start=1;
k_end=25;
% Maximum observed frequency modulation
psi_thresh=100/Fs*2*pi;
load([path,'.dat']);
for n=k_start:k_end
    X_a=a(:,n);
    % only plot maxima under nyquist
    [idmx,Y]=ppvp(X_a(1:(end/2)),0.9,'max',0.50,2);
%    d=diff(X_a);
    % only plot maxima greater than threshold
    idmx=idmx(find(X_a(idmx)>thresh));
    % only plot maxima with a reasonable frequency slope
    idmx=idmx(find(abs(psi(idmx,n))<psi_thresh));
    X_max=X_a(idmx);
    w_max=w(idmx,n);
    n_max=idmx;
    t_max=t(idmx,n);
    psi_max=psi(idmx,n);
    t0=-H/2*ones(length(idmx),1);
    t1=H/2*ones(length(idmx),1);
    w0=w_max+psi_max.*t0;
    w1=w_max+psi_max.*t1;
%    h=scatter((n-1)*ones(length(idmx),1)*H+t_max,w_max/(2*pi)*Fs,2,...
%        clip(X_max*ones(1,3)));
    % Note that we do not use the reassigned time
    h=scatter((n-1)*ones(length(idmx),1)*H,w_max/(2*pi)*Fs,2,...
        clip(X_max*ones(1,3)));
    set(h,'linewidth',2);
%    plot(([t0+t_max,t1+t_max]+(n-1)*ones(length(idmx),2)*H)',[w0,w1]'/(2*pi)*Fs);
    plot(([t0,t1]+(n-1)*ones(length(idmx),2)*H)',[w0,w1]'/(2*pi)*Fs,'b','linewidth',2);
    % compare with theoretical value of slope
    f_th=sp.f0+sp.A_fm*sin(2*pi*sp.f_fm*((n-1)*H/Fs)+sp.phi_fm)...
        +2*pi*sp.A_fm*sp.f_fm*cos(2*pi*sp.f_fm*(n-1)*H/Fs+sp.phi_fm)...
        *([t0(1) t1(1)])/Fs;
    f_th=(sp.k_B)'*f_th;
%    f_th=sp.f0+sp.A_fm*sin(2*pi*sp.f_fm*((n-1)*H/Fs)+sp.phi_fm);
%    f_th=[f_th f_th];
    plot(([t0(1),t1(1)].+(n-1)*ones(size(f_th))*H)',f_th','r','linewidth',2);
end
hold off;
