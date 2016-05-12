% Compute theoretical values with harm_sines_rp
% Try to separate them using GMM
% Here spurious data are added to see how they distrupt the analysis.
% The resulting partials (minus the spurious ones) are saved so that the
% partials can be grouped accross time frames by a different algorithm.
% The file produced by this is not compatible with those produced by
% hsrp_test_8.m
clear;
fname=sprintf('/tmp/hsrpc_%s.mat',datestr(now(),'yyyymmddTHHMMSSFFF'));
printf('%s\n',fname)
% Number of iterations of GMM algorithm
B_gmm=50;
% Number of histogram bins
N_b=40;
% Histrogram bin count threshold
H_th=7;
N_pxm=2;
Pxm=cell(N_pxm);
out=cell(3);
[Pxm{1},opt]=harm_sines_rp(struct(
    'K',10,
    'T',.5,
    'H',256,
    'f0',440*2^((60-69)/12),
    'w_no',0.001,
    'T_60',0.5,
    'mu_no',0.001,
    'f_fm',3,
    'psi_no',0.001));
out{1,3}=opt;
[Pxm{2},opt]=harm_sines_rp(struct(
    'K',10,
    'T',.5,
    'H',256,
    'f0',440*2^((61-69)/12),
    'w_no',0.001,
    'T_60',0.75,
    'mu_no',0.001,
    'phi_fm',3,
    'f_fm',2,
    'psi_no',0.001));
out{2,3}=opt;
% percentage of spurious peaks added in relation to number of real peaks.
spur_no=0.5;
% range of fake w parameters
w_no_rng=[0 pi/2];
% range of fake psi parameters
psi_no_rng=[-1e-4 1e-4];
% range of fake mu parameters
mu_no_rng=[-1e-4 1e-4];
phi_no_rng=[-pi pi];
A_no_rng=[1e-4 1];
for n=1:length(Pxm{1})
    X=[]; % storage for random variables
    X_plt=[]; % Storage for plotting information
    m=(n-1)*opt.H;
    for n_pxm=1:N_pxm
        sp=Pxm{n_pxm}{n};
        n0=-opt.H/2;
        n1=opt.H/2;
        X=[X;[(sp.psi_r./sp.w_r)(:),sp.mu_r(:)]];
        % Vector of values to be plotted.
        X_plt=[X_plt;[sp.w_r,sp.psi_r,abs(sp.X_r_),arg(sp.X_r_),sp.mu_r]];
    end
    % Store the true values (uncorrupted)
    out{n,4}=X_plt;
    N_fake=round(size(X,1)*spur_no);
    w_fake=unifrnd(w_no_rng(1),w_no_rng(2),N_fake,1);
    psi_fake=unifrnd(psi_no_rng(1),psi_no_rng(2),N_fake,1);
    mu_fake=unifrnd(mu_no_rng(1),mu_no_rng(2),N_fake,1);
    phi_fake=unifrnd(phi_no_rng(1),phi_no_rng(2),N_fake,1);
    A_fake=unifrnd(A_no_rng(1),A_no_rng(2),N_fake,1);
    X_fake=A_fake.*exp(j.*phi_fake);
    X=[X;[psi_fake./w_fake,mu_fake]];
    X_plt=[X_plt;[w_fake,psi_fake,A_fake,phi_fake,mu_fake]];
    % Permute the data vectors to prove algorithm works.
    perm_idx=randperm(size(X,1))(:);
    X=X(perm_idx,:);
    X_plt=X_plt(perm_idx,:);
    % First PC of these as realizations of 2-dimensional random variable
    [A_pca,l_pca]=pca_ne(X,'corr');
    [Y,Xdi,Xki]=hist_thresh(A_pca(:,1),N_b,H_th);
    % First guess for means and standard deviations
    s0=zeros(1,1,N_pxm);
    [m0,s0(1,1,:),b_avg,P_modes]=hist_sect(Y,N_pxm,N_b);
    % First guess for weights
    w0=ones(N_pxm,1)/N_pxm;
    [C_,P_,mu_,S_,W_]=gmm(Y,m0,(s0/3).^2,w0,B_gmm);
    out{n,1}=X;
    out{n,2}=X_plt;
    out{n,5}=A_pca;
    % Indices of data to discard
    out{n,6}=Xdi;
    % Indices of data to keep
    out{n,7}=Xki;
    % Classifications of kept data
    out{n,8}=C_;
end
save('-6',fname,'out')
