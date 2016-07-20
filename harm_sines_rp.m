function [Pxm,opt]=harm_sines_rp(opt=struct())
% Directly compute the theoretical reassigned parameters
% default synthesis parameters
sp=struct();
% Sampling rate in Hz
sp.Fs=16000;
% Length of signal in seconds
sp.T=1.5;
% Fundamental
sp.f0=440*2^((67-69)/12);
% Number of harmonics
sp.K=10;
% Inharmonicity coefficient
%sp.B=exp(2.54*log(sp.f0)-24.6); % the just noticeable inharmonic coefficient
sp.B=0.001; % an ad hoc value
% Amplitude based on harmonic number
sp.A_k_60=20; % The harmonic number that is 60dB lower than the first
% Initial phase
sp.phi=0;
% Initial FM phase
sp.phi_fm=0;
% Amplitude coefficient of FM
sp.A_fm=sp.f0*2^(0.5/12)-sp.f0;
% Frequency coefficient of FM
sp.f_fm=2;
% Time in seconds until amplitude of partial has dropped by 60dB
sp.T_60=.75;
% Time in seconds to when attack function reaches maximum
sp.T_max=0.1;
% Attack method
sp.A_method='FOF';
% Theoretical size of analysis frame (actual size will sp.L + 1)
sp.L=1024;
% Theoretical hop size
sp.H=256;
% For the noise amounts, this value can be interpreted as: 68% of the time, the
% value will be within this percentage of its true value, where 0.01 is 1% etc.
% frequency noise variance
sp.w_no=0;
% Same as for frequency noise
% amplitude noise variance
sp.A_no=0;
% For the frequency slope noise variance, this is a percentage of the minimum
% and maximum theoretical slope of the fundamental harmonic. No more than this
% value will be added or subtracted from the theoretical psi value 68% of the
% time
% frequency slope noise variance
sp.psi_no=0;
% For the amplitude slope noise variance, this is the amount in dB that the
% amplitude after 1 hop's worth of samples will deviate from the theoretical
% amplitude, had there not been any noise
% amplitude slope noise variance
sp.mu_no=0;

% Check what fields are present and replace with defaults
for fn=fieldnames(sp)'
    fn=char(fn);
    if(~isfield(opt,fn))
        opt.(fn)=sp.(fn);
    end
end

% Length of signal in samples
opt.N=opt.Fs*opt.T;
% Harmonic numbers
k=1:opt.K;
% Harmonic numbers adjusted for inharmonicity
opt.k_B=k.*sqrt(1+opt.B*k.^2);
%a_k_60=log(10^(-3))/log(10)/opt.A_k_60;
a_k_60=log(10^(-3))/opt.A_k_60;
% use original harmonic numbers, divide by 2 because
% we store only one half of cosine in spectrum (not its conjugate phasor)
opt.A_k=exp(a_k_60*k);
% amplitude coefficient
opt.a_60=log(10^(-3))/opt.T_60;
% psi stddev calculation
psi_sd=2*(2*pi)^2*opt.A_fm*opt.f_fm/(opt.Fs^2)*opt.psi_no;
% mu stddev calculation
mu_sd=opt.mu_no*log(10)/(20*opt.H);
% sample indices
n=(0:opt.H:(opt.N-1));
n=n(:);
rp=struct();
t=n/opt.Fs;
dt=t(2)-t(1);
correlated_measurements=0;
Pxm=cell();
n=1;
for t_=t'
    switch opt.A_method
    case 'FOF'
        % using FOF functions
        if (t_ < opt.T_max)
            % For attack portion we use least-squares fit to estimate mu (AM) parameter
            n_=[-opt.L/2:opt.L/2]';
            t__=t_+n_/opt.Fs;
            t__=t__(:);
            s_=exp(opt.a_60*t__)*opt.A_k.*cos(t__/opt.T_max*pi/2-pi/2).^2;
            % mu (AM) parameters are in row 2, log of initial amplitude in row 1
            % The regressors (?) are with repect to sample numbers rather than
            % time because this is how the parameters are computed by RM
            th=ols(log(s_),[ones(length(t__),1),n_]);
            rp.mu_r=th(2,:)';
            if correlated_measurements==1
                error('Currently unknown how to implement.')
%            % Add noise
%%            murand=randn(length(rp.mu_r),1).*(mu_sd);
%            murand=randn(length(rp.mu_r),1).*(opt.mu_no)*0.;
%            rp.mu_r.+=murand;
%            alphrand=murand*dt+randn(length(th(1,:)),1)*opt.A_no;
%            alph_r=th(1,:)'+alphrand;
%%            rp.X_r_=exp(th(1,:)').*(1+randn(length(th(1,:)),1)*opt.A_no);
%            rp.X_r_=exp(alph_r);
            else
                murand=randn(length(rp.mu_r),1).*(opt.mu_no);
                rp.mu_r.+=murand;
                alphrand=randn(length(th(1,:)),1)*opt.A_no;
                alph_r=th(1,:)'+alphrand;
                rp.X_r_=exp(alph_r);
            end

        else
            % Divided by Fs because opt.a_60 is w.r.t. the time in seconds not the
            % sample number
            rp.mu_r=ones(opt.K,1)*opt.a_60/opt.Fs;
            % Add noise
%            murand=randn(length(rp.mu_r),1).*(mu_sd);
            if correlated_measurements==1
                error('Currently unknown how to implement.')
%            murand=randn(length(rp.mu_r),1).*(opt.mu_no)*0.;
%            rp.mu_r.+=murand;
%            alphrand=murand*dt+randn(opt.K,1)*opt.A_no;
%            alph_r=opt.a_60*t_+alphrand+log(opt.A_k');
%            rp.X_r_=exp(alph_r);
%            % Add noise
%            % The value where the AM parameter is multiplied by 0 is at the
%            % centre of the window, which is then just the amplitude as if there
%            % were no AM
%            %rp.X_r_=exp(opt.a_60*t_)*(opt.A_k'.*(1+randn(opt.K,1)*opt.A_no));
            else
                murand=randn(length(rp.mu_r),1).*(opt.mu_no);
                rp.mu_r.+=murand;
                alphrand=randn(opt.K,1)*opt.A_no;
                alph_r=opt.a_60*t_+alphrand+log(opt.A_k');
                rp.X_r_=exp(alph_r);
            end

        end
    otherwise
        error(sprintf('Bad A_method %s\n',opt.A_method));
    end
    % Reassigned frequency parameters
    if correlated_measurements==1
        psirand=randn(opt.K,1)*psi_sd;
        wrand=psirand*dt+randn(opt.K,1)*opt.w_no;
        phirand=wrand*dt;
        rp.w_r=2*pi*(opt.f0+opt.A_fm*sin(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B'/opt.Fs;
        rp.w_r+=wrand;
        rp.psi_r=(2*pi)^2*opt.A_fm*opt.f_fm*cos(2*pi*opt.f_fm*t_+opt.phi_fm)*opt.k_B'/(opt.Fs^2);
        rp.psi_r.+=psirand;
        % Initial phase
        rp.phi_r=opt.phi+2*pi*(opt.f0*t_-opt.A_fm/(2*pi*opt.f_fm)*cos(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B';
        rp.phi_r+=phirand;
        rp.X_r_.*=exp(j*rp.phi_r);
    else
        psirand=randn(opt.K,1)*psi_sd;
        wrand=randn(opt.K,1)*opt.w_no;
        phirand=0;
        rp.w_r=2*pi*(opt.f0+opt.A_fm*sin(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B'/opt.Fs;
        rp.w_r+=wrand;
        rp.psi_r=(2*pi)^2*opt.A_fm*opt.f_fm*cos(2*pi*opt.f_fm*t_+opt.phi_fm)*opt.k_B'/(opt.Fs^2);
        rp.psi_r.+=psirand;
        % Initial phase
        rp.phi_r=opt.phi+2*pi*(opt.f0*t_-opt.A_fm/(2*pi*opt.f_fm)*cos(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B';
        rp.phi_r+=phirand;
        rp.X_r_.*=exp(j*rp.phi_r);
    end
    Pxm{n}=rp;
    n=n+1;
end
