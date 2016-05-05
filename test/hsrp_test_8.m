% Compute theoretical values with harm_sines_rp 
% Add spurious data
% combine parameters in each frame
% save result (for loading into other programs)
clear;
fname=sprintf('/tmp/hsrp_%s.mat',datestr(now(),'yyyymmddTHHMMSSFFF'));
N_pxm=2;
Pxm=cell(N_pxm);
out=cell(3);
[Pxm{1},opt]=harm_sines_rp(struct(
    'K',5,
    'T',0.125,
    'H',256,
    'f0',440*2^((60-69)/12),
    'w_no',0.001,
    'T_60',0.5,
    'mu_no',0.001,
    'f_fm',3,
    'psi_no',0.001));
out{1,3}=opt;
[Pxm{2},opt]=harm_sines_rp(struct(
    'K',5,
    'T',0.125,
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
spur_no=2;
% range of fake w parameters
w_no_rng=[0 pi/4];
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
    out{n,1}=X;
    out{n,2}=X_plt;
end
save('-6',fname,'out')
