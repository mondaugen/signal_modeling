% Compute theoretical values with harm_sines_rp and then plot them
[Pxm,opt]=harm_sines_rp();
newplot(figure(1));
figure(1);
hold on;
for n=1:length(Pxm)
    sp=Pxm{n};
    h=scatter(n*ones(length(sp.w_r),1),sp.w_r,[],'k');
    set(h,'linewidth',1);
    n0=-opt.H/2;
    n1=opt.H/2;
    w0=sp.psi_r*n0;
    w1=sp.psi_r*n1;
    plot([-0.5,0.5]+n,[w0,w1]+sp.w_r,'k');
end
hold off;

