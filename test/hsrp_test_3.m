% Compute theoretical values with harm_sines_rp and then plot the frequency
% modulation / frequency and amplitude modulation / amplitude
clear;
N_pxm=2;
Pxm=cell(N_pxm);
[Pxm{1},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((60-69)/12),
    'w_no',0.1,
    'T_60',0.5,
    'mu_no',0.1,
    'f_fm',3,
    'psi_no',0.1));
[Pxm{2},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((61-69)/12),
    'w_no',0.1,
    'T_60',0.75,
    'mu_no',0.1,
    'phi_fm',.8,
    'f_fm',2,
    'psi_no',0.1));
f1=newplot(figure(1));
f2=newplot(figure(2));
f3=newplot(figure(3));
f4=newplot(figure(4));
hold(f1,'on');
hold(f2,'on');
hold(f3,'on');
hold(f4,'on');
A_pca=cell();
for n=1:length(Pxm{1})
    X=[]; % storage for random variables
    for n_pxm=1:N_pxm
        sp=Pxm{n_pxm}{n};
        m=(n-1)*opt.H;
        figure(1);
        % Frequency modulation / frequency
        h=scatter(m*ones(length(sp.w_r),1),sp.psi_r./sp.w_r,[],'k');
        set(h,'linewidth',1);
        n0=-opt.H/2;
        n1=opt.H/2;
        figure(2);
        % Amplitude modulation
        h=scatter(m*ones(length(sp.mu_r),1),sp.mu_r,[],'k');
        set(h,'linewidth',1);
        X=[X;[(sp.psi_r./sp.w_r)(:),sp.mu_r(:)]];
    end
    % First PC of these as realizations of 2-dimensional random variable
    [A_pca{n},l_pca]=pca_ne(X,'corr');
    figure(3);
    h=scatter(m*ones(length(A_pca{n}(:,1)),1),A_pca{n}(:,1),[],'k');
    set(h,'linewidth',1);
    figure(4);
    % Second PC of these as realizations of 2-dimensional random variable
    h=scatter(m*ones(length(A_pca{n}(:,2)),1),A_pca{n}(:,2),[],'k');
    set(h,'linewidth',1);
end
%for n=1:length(Pxm2)
%    sp=Pxm2{n};
%    m=(n-1)*opt.H;
%    figure(1);
%    h=scatter(m*ones(length(sp.w_r),1),sp.psi_r./sp.w_r,[],'k');
%    set(h,'linewidth',1);
%    n0=-opt.H/2;
%    n1=opt.H/2;
%    w0=sp.psi_r*n0;
%    w1=sp.psi_r*n1;
%%    plot([-0.5,0.5]*opt.H+m,[w0,w1]+sp.w_r,'k');
%    figure(2);
%%    h=scatter(m*ones(length(sp.mu_r),1),sp.mu_r./log(abs(sp.X_r_)),[],'k');
%    h=scatter(m*ones(length(sp.mu_r),1),sp.mu_r,[],'k');
%    set(h,'linewidth',1);
%    a0=sp.mu_r*n0;
%    a1=sp.mu_r*n1;
%%    plot([-0.5,0.5]*opt.H+m,[a0,a1]+log(abs(sp.X_r_)),'k');
%end
%M=opt.H*length(Pxm);
%m=(0:(M-1));
%t_=m/opt.Fs;
%t_=t_(:);
%w_r=2*pi*(opt.f0+opt.A_fm*sin(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B/opt.Fs;
%figure(1);
%plot(m,w_r,'--');
%a_r=log(exp(opt.a_60*t_)*opt.A_k.*(cos(t_/opt.T_max*pi/2-pi/2).^2.*(t_<=opt.T_max)+(t_>opt.T_max)));
%figure(2);
%plot(m,a_r,'--');
hold(f1,'off');
hold(f2,'off');
hold(f3,'off');
hold(f4,'off');
