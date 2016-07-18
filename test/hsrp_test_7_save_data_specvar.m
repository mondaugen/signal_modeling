% Compute theoretical values with harm_sines_rp and then plot the frequency
% modulation / frequency and amplitude modulation / amplitude
% Here we undertake the additional step of trying to separate them using GMM
% We then plot the separated source components on a time-frequency plot with
% frequency slope vector fields.
% Here spurious data are added to see how they distrupt the analysis.
% No plotting, save data for plotting with matplotlib.
% Specify parameter noise variance as argument to script
clear -x _int_args;
if (exist('_int_args'))
    argv_=_int_args;
else
    argv_=argv();
end
if length(argv_) != 2
    error('Specify desired variance')
else
    no_f=str2num(argv_{1});
    no_a=str2num(argv_{2});
end
    
datoutpath=[getenv('HOME'),'/Documents/development/masters_thesis/reports/plots/'];
datoutpath=[datoutpath,'hsrp_test_7_',argv_{1},'.dat'];
datoutpath
% Colours for plotting different categories
clrs={"black","blue","cyan","green","magenta","red","yellow"};
% Number of iterations of GMM algorithm
B_gmm=50;
N_pxm=2;
Pxm=cell(N_pxm);
[Pxm{1},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((60-69)/12),
    'w_no',no_f,
    'T_60',0.5,
    'mu_no',no_a,
    'A_no',no_a,
    'f_fm',3,
    'psi_no',no_f));
[Pxm{2},opt]=harm_sines_rp(struct(
    'K',20,
    'T',0.5,
    'H',256,
    'f0',440*2^((61-69)/12),
    'w_no',no_f,
    'T_60',0.75,
    'mu_no',no_a,
    'A_no',no_a,
    'phi_fm',.8,
    'f_fm',2,
    'psi_no',no_f));
% percentage of spurious peaks added in relation to number of real peaks.
spur_no=0.25;
% range of fake w parameters
w_no_rng=[0 pi];
% range of fake psi parameters
psi_no_rng=[-1e-4 1e-4];
% range of fake mu parameters
mu_no_rng=[-1e-3 1e-3];
% Number of histogram bins
N_b=20;
% Histrogram bin count threshold
H_th=7;
A_pca=cell();
datout=cell();
for n=1:length(Pxm{1})
    datout{n}=struct();
    X=[]; % storage for random variables
    X_plt=[]; % Storage for plotting information
    lA_plt=[]; % Storage for log amplitude
    m=(n-1)*opt.H;
    for n_pxm=1:N_pxm
        sp=Pxm{n_pxm}{n};
        n0=-opt.H/2;
        n1=opt.H/2;
        X=[X;[(sp.psi_r./sp.w_r)(:),sp.mu_r(:)]];
        % Vector of values to be plotted.
        X_plt=[X_plt;[sp.w_r,sp.psi_r]];
        lA_plt=[lA_plt;[log(abs(sp.X_r_))]];
    end
    datout{n}.X_orig=X;
    datout{n}.X_plt_orig=X_plt;
    datout{n}.lA_plt_orig=lA_plt;
    N_fake=round(size(X,1)*spur_no);
    w_fake=unifrnd(w_no_rng(1),w_no_rng(2),N_fake,1);
    psi_fake=unifrnd(psi_no_rng(1),psi_no_rng(2),N_fake,1);
    mu_fake=unifrnd(mu_no_rng(1),mu_no_rng(2),N_fake,1);
    % Amplitude is not considered in classification so just put some junk
    lA_plt_fake=unifrnd(-80,0,N_fake,1);
    X=[X;[psi_fake./w_fake,mu_fake]];
    X_plt=[X_plt;[w_fake,psi_fake]];
    lA_plt=[lA_plt;lA_plt_fake];
    datout{n}.X_fake=[psi_fake./w_fake,mu_fake];
    datout{n}.X_plt_fake=[w_fake,psi_fake];
    datout{n}.lA_plt_fake=lA_plt_fake;
    % Frequency modulation / frequency
%    h=scatter(f1,m*ones(size(X,1),1),X(:,1),[],'k');
%    set(h,'linewidth',1);
    % Amplitude modulation
%    h=scatter(f2,m*ones(size(X,1),1),X(:,2),[],'k');
%    set(h,'linewidth',1);
    % Permute the data vectors to prove algorithm works.
    perm_idx=randperm(size(X,1))(:);
    X=X(perm_idx,:);
    X_plt=X_plt(perm_idx,:);
    lA_plt=lA_plt(perm_idx,:);
    datout{n}.X=X;
    datout{n}.X_plt=X_plt;
    datout{n}.lA_plt=lA_plt;
    % First PC of these as realizations of 2-dimensional random variable
    [A_pca{n},l_pca]=pca_ne(X,'corr');
    datout{n}.A_pca=A_pca{n};
    [Y,Xdi,Xki]=hist_thresh(A_pca{n}(:,1),N_b,H_th);
    % First guess for means and standard deviations
    s0=zeros(1,1,N_pxm);
    [m0,s0(1,1,:),b_avg,P_modes]=hist_sect(Y,N_pxm,N_b);
    datout{n}.m0=m0;
    datout{n}.b_avg=b_avg;
    datout{n}.s0=(s0/3).^2;
    % First guess for weights
    w0=ones(N_pxm,1)/N_pxm;
    [C_,P_,mu_,S_,W_]=gmm(Y,m0,(s0/3).^2,w0,B_gmm);
%   Instead of Gaussian mixture models, we just pick points within the
%   discovered boundaries as belonging to a given source
%    C_=zeros(size(A_pca{n}(:,1)));
%    C_(find(A_pca{n}(:,1) < b_avg(2)))=1;
%    C_(find(A_pca{n}(:,1) >= b_avg(2)))=2;
    for c=1:N_pxm
        ci_=find(C_==c);
        %h=scatter(f3,m*ones(length(Y(ci_,1)),1),Y(ci_,1),[],clrs{c+1});
        %set(h,'linewidth',1);
    end
    % plot uncategorized values
%    h=scatter(f3,m*ones(length(A_pca{n}(Xdi,1)),1),A_pca{n}(Xdi,1),[],clrs{1});
%    set(h,'linewidth',1);
    % plot first guess for means
%    h=scatter(f3,m*ones(length(m0),1),m0,[],'g');
%    set(h,'linewidth',1);
    % plot area boundaries
%    h=scatter(f3,m*ones(length(b_avg),1),b_avg,[],'r');
%    set(h,'linewidth',1);
%    title(f3,'First principal component (PC)');
%    ylabel(f3,'PC value');
%    xlabel(f3,'Sample number');
    w0=X_plt(:,1)+X_plt(:,2)*opt.H*(-0.5);
    w1=X_plt(:,1)+X_plt(:,2)*opt.H*(0.5);
    datout{n}.w0=w0;
    datout{n}.w1=w1;
    clr_=zeros(size(X_plt,1),1);
    clr_(Xki(find(C_==1)))='b';
    clr_(Xki(find(C_==2)))='c';
    clr_(Xdi)='k';
    clr_=char(clr_);
    datout{n}.clr_=clr_;
    %h=plot(f4,m+opt.H*[-0.5 0.5],[w0,w1],clr_);
%    set(h,'linewidth',2);
%    title('Partial groupings');
%    xlabel('Sample number');
%    ylabel('Frequency (rad/s)');
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
save('-mat',datoutpath,'datout')
%hold(f1,'off');
%hold(f2,'off');
%hold(f3,'off');
%hold(f4,'off');
%set(f1,'visible','on');
%set(f2,'visible','on');
%set(f3,'visible','on');
%set(f4,'visible','on');
