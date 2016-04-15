% Compute theoretical values with harm_sines_rp and then plot them
[Pxm,opt]=harm_sines_rp(struct('T',2,'H',256));
newplot(figure(1));
newplot(figure(2));
hold on;
for n=1:length(Pxm)
    sp=Pxm{n};
    m=(n-1)*opt.H;
    figure(1);
    h=scatter(m*ones(length(sp.w_r),1),sp.w_r,[],'k');
    set(h,'linewidth',1);
    n0=-opt.H/2;
    n1=opt.H/2;
    w0=sp.psi_r*n0;
    w1=sp.psi_r*n1;
    plot([-0.5,0.5]*opt.H+m,[w0,w1]+sp.w_r,'k');
    figure(2);
    h=scatter(m*ones(length(sp.mu_r),1),log(abs(sp.X_r_)),[],'k');
    set(h,'linewidth',1);
    a0=sp.mu_r*n0;
    a1=sp.mu_r*n1;
    plot([-0.5,0.5]*opt.H+m,[a0,a1]+log(abs(sp.X_r_)),'k');
end
M=opt.H*length(Pxm);
m=(0:(M-1));
t_=m/opt.Fs;
t_=t_(:);
w_r=2*pi*(opt.f0+opt.A_fm*sin(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B/opt.Fs;
figure(1);
plot(m,w_r,'--');
a_r=log(exp(opt.a_60*t_)*opt.A_k.*(cos(t_/opt.T_max*pi/2-pi/2).^2.*(t_<=opt.T_max)+(t_>opt.T_max)));
figure(2);
plot(m,a_r,'--');
hold off;
