% Compute theoretical values with harm_sines_rp and then plot the frequency
% modulation / frequency and amplitude modulation / amplitude
% Here we undertake the additional step of trying to separate them using GMM
% We then plot the separated source components on a time-frequency plot with
% frequency slope vector fields.
% Here little noise is added to show clearly the method works. It works for more
% noise, but the sources in the time-frequency plot are unrecognizable.
clear;
% Colours for plotting different categories
clrs={"blue", "black", "cyan", "green", "magenta", "red","yellow"};
% Number of iterations of GMM algorithm
B_gmm=50;
N_pxm=2;
Pxm=cell(N_pxm);
[Pxm{1},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((60-69)/12),
    'w_no',0.001,
    'T_60',0.5,
    'mu_no',0.001,
    'f_fm',3,
    'psi_no',0.001));
[Pxm{2},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((61-69)/12),
    'w_no',0.001,
    'T_60',0.75,
    'mu_no',0.001,
    'phi_fm',.8,
    'f_fm',2,
    'psi_no',0.001));
f1=newplot(figure(1));
set(f1,'visible','off');
f2=newplot(figure(2));
set(f2,'visible','off');
f3=newplot(figure(3));
set(f3,'visible','off');
f4=newplot(figure(4));
set(f4,'visible','off');
hold(f1,'on');
hold(f2,'on');
hold(f3,'on');
hold(f4,'on');
A_pca=cell();
for n=1:length(Pxm{1})
    X=[]; % storage for random variables
    X_plt=[]; % Storage for plotting information
    m=(n-1)*opt.H;
    for n_pxm=1:N_pxm
        sp=Pxm{n_pxm}{n};
        % Frequency modulation / frequency
        h=scatter(f1,m*ones(length(sp.w_r),1),sp.psi_r./sp.w_r,[],'k');
        set(h,'linewidth',1);
        n0=-opt.H/2;
        n1=opt.H/2;
        % Amplitude modulation
        h=scatter(f2,m*ones(length(sp.mu_r),1),sp.mu_r,[],'k');
        set(h,'linewidth',1);
        X=[X;[(sp.psi_r./sp.w_r)(:),sp.mu_r(:)]];
        % Vector of values to be plotted.
        X_plt=[X_plt;[sp.w_r,sp.psi_r]];
    end
    % Permute the data vectors to prove algorithm works.
    perm_idx=randperm(size(X,1))(:);
    X=X(perm_idx,:);
    X_plt=X_plt(perm_idx,:);
    % First PC of these as realizations of 2-dimensional random variable
    [A_pca{n},l_pca]=pca_ne(X,'corr');
    % First guess for means and standard deviations
    s0=zeros(1,1,N_pxm);
    [m0,s0(1,1,:)]=hist_sect(A_pca{n}(:,1),N_pxm);
    %[m0,s0(1,1,:),b_avg,P_modes]=hist_sect(A_pca{n}(:,1),N_pxm);
    %% Instead use local maxima in histogram as initial means
    %m0=P_modes;
    % First guess for weights
    w0=ones(N_pxm,1)/N_pxm;
    [C_,P_,mu_,S_,W_]=gmm(A_pca{n}(:,1),m0,s0*0.001,w0,B_gmm);
    for c=1:N_pxm
        ci_=find(C_==c);
        h=scatter(f3,m*ones(length(A_pca{n}(ci_,1)),1),A_pca{n}(ci_,1),[],clrs{c});
        set(h,'linewidth',1);
    end
    % plot first guess for means
    h=scatter(f3,m*ones(length(m0),1),m0,[],'g');
    set(h,'linewidth',1);
    w0=X_plt(:,1)+X_plt(:,2)*opt.H*(-0.5);
    w1=X_plt(:,1)+X_plt(:,2)*opt.H*(0.5);
    clr_=zeros(size(X_plt,1),1);
    clr_(find(C_==1))='k';
    clr_(find(C_==2))='b';
    clr_=char(clr_);
    plot(f4,m+opt.H*[-0.5 0.5],[w0,w1],clr_);
    %% First guess for means and standard deviations
    %[m0,s0(1,1,:)]=hist_sect(A_pca{n}(:,2),N_pxm);
    %% First guess for weights
    %w0=ones(N_pxm,1)/N_pxm;
    %[C_,P_,mu_,S_,W_]=gmm(A_pca{n}(:,1),m0,s0,w0,B_gmm);
    %for c=1:N_pxm
    %    ci_=find(C_==c);
    %    h=scatter(f4,m*ones(length(A_pca{n}(ci_,2)),1),A_pca{n}(ci_,2),[],clrs{c});
    %    set(h,'linewidth',1);
    %end
    %% plot first guess for means
    %h=scatter(f4,m*ones(length(m0),1),m0,[],'g');
    %set(h,'linewidth',1);
end
hold(f1,'off');
hold(f2,'off');
hold(f3,'off');
hold(f4,'off');
set(f1,'visible','on');
set(f2,'visible','on');
set(f3,'visible','on');
set(f4,'visible','on');
