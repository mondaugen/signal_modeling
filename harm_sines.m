function [x,opt]=harm_sines(opt=struct())
%synthesis parameters
% Sampling rate in Hz
sp.Fs=16000;
% Length of signal
sp.N=sp.Fs*1.5;
% Fundamental
sp.f0=440*2^((67-69)/12);
% Number of harmonics
sp.K=10;
% Inharmonicity coefficient
%sp.B=exp(2.54*log(sp.f0)-24.6); % For now, the just noticeable inharmonic coefficient
sp.B=0.001;
% Amplitude based on harmonic number
sp.A_k_60=20; % The harmonic number that is 60dB lower than the first
% Initial phase
sp.phi=0;
% Initial FM phase
sp.phi_fm=0;
% Time in seconds until amplitude of partial has dropped by 60dB
%sp.T_60=1;
sp.T_60=.75;
% Time in seconds to when attack function reaches maximum
sp.T_max=0.1;
sp.A_method='FOF';

% Check what fields are present and replace with defaults
for fn=fieldnames(sp)'
    fn=char(fn);
    if(~isfield(opt,fn))
        opt.(fn)=sp.(fn);
    end
end

% amplitude coefficient
a_60=log(10^(-3))/log(10)/opt.T_60;
% Harmonic numbers
k=1:opt.K;
% Harmonic numbers adjusted for inharmonicity
opt.k_B=k.*sqrt(1+opt.B*k.^2);
a_k_60=log(10^(-3))/log(10)/opt.A_k_60;
opt.A_k=exp(a_k_60*k); % use original harmonic numbers
% sample indices
n=(0:(opt.N-1));
n=n(:);
t=n/opt.Fs;

switch opt.A_method
case 'FOF'
    % using FOF functions
    A=exp(a_60*t)*opt.A_k.*(cos(t/opt.T_max*pi/2-pi/2).^2.*(t<=opt.T_max)+(t>opt.T_max));
case 'AR1'
    % Using first order AR model
    C_max=0.99;
    n_max=opt.Fs*opt.T_max;
    a_T_max=-(1-C_max)^(1/(n_max+1));
    A=filter([1+a_T_max],[1 a_T_max],exp(a_60*t)*opt.A_k,[],1);
otherwise
    error('Bad A_method');
end
% Amplitude coefficient of FM
opt.A_fm=opt.f0*2^(1/12)-opt.f0;
% Frequency coefficient of FM
opt.f_fm=2;
x=A.*cos(2*pi*(opt.f0*t-opt.A_fm/(2*pi*opt.f_fm)*cos(2*pi*opt.f_fm*t+opt.phi_fm))*opt.k_B+opt.phi);
x=sum(x,2);
