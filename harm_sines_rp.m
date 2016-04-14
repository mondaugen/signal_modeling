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
sp.A_fm=sp.f0*2^(2/12)-sp.f0;
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
a_k_60=log(10^(-3))/log(10)/opt.A_k_60;
% use original harmonic numbers, divide by 2 because
% we store only one half of cosine in spectrum (not its conjugate phasor)
opt.A_k=exp(a_k_60*k)*0.5;
% amplitude coefficient
a_60=log(10^(-3))/log(10)/opt.T_60;
% sample indices
n=(0:opt.H:(opt.N-1));
n=n(:);
rp=struct();
t=n/opt.Fs;
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
               s_=exp(a_60*t__)*opt.A_k.*cos(t__/opt.T_max*pi/2-pi/2).^2;
               % mu (AM) parameters are in row 2, log of initial amplitude in row 1
               % The regressors (?) are with repect to sample numbers rather than
               % time because this is how the parameters are computed by RM
               th=ols(log(s_),[ones(length(t__),1),n_]);
               rp.mu_r=th(2,:)';
               rp.X_r_=exp(th(1,:)');
           else
               % Divided by Fs because a_60 is w.r.t. the time in seconds not the
               % sample number
               rp.mu_r=ones(opt.K,1)*a_60/opt.Fs;
               % The value where the AM parameter is multiplied by 0 is at the
               % centre of the window, which is then just the amplitude as if there
               % were no AM
               rp.X_r_=exp(a_60*t_)*opt.A_k';
           end
    otherwise
        error(sprintf('Bad A_method %s\n',opt.A_method));
    end
    % Reassigned frequency parameters
    rp.w_r=(2*pi*opt.f0+opt.A_fm*sin(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B'/opt.Fs;
    rp.psi_r=(2*pi)^2*opt.A_fm*opt.f_fm*cos(2*pi*opt.f_fm*t_+opt.phi_fm)*opt.k_B'/(opt.Fs^2);
    % Initial phase
    rp.phi_r=opt.phi+2*pi*(opt.f0*t_-opt.A_fm/(2*pi*opt.f_fm)*cos(2*pi*opt.f_fm*t_+opt.phi_fm))*opt.k_B';
    rp.X_r_.*=exp(j*rp.phi_r);
    Pxm{n}=rp;
    n=n+1;
end
