function [x,opt]=harm_sines(opt=struct())
%synthesis parameters
% Sampling rate in Hz
sp.Fs=16000;
% Length of signal
sp.N=sp.Fs*1.5;
% sample indices
n=(0:(sp.N-1));
n=n(:);
t=n/sp.Fs;
% Fundamental
sp.f0=440*2^((67-69)/12);
% Number of harmonics
sp.K=10;
% Inharmonicity coefficient
%sp.B=exp(2.54*log(sp.f0)-24.6); % For now, the just noticeable inharmonic coefficient
sp.B=0.001;
% Harmonic numbers
k=1:sp.K;
% Harmonic numbers adjusted for inharmonicity
sp.k_B=k.*sqrt(1+sp.B*k.^2);
% Amplitude based on harmonic number
sp.A_k_60=20; % The harmonic number that is 60dB lower than the first
a_k_60=log(10^(-3))/log(10)/sp.A_k_60;
sp.A_k=exp(a_k_60*k); % use original harmonic numbers
% Initial phase
sp.phi=0;
% Initial FM phase
sp.phi_fm=0;
% Time in seconds until amplitude of partial has dropped by 60dB
%sp.T_60=1;
sp.T_60=.75;
% amplitude coefficient
a_60=log(10^(-3))/log(10)/sp.T_60;
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

switch sp.A_method
case 'FOF'
    % using FOF functions
    A=exp(a_60*t)*sp.A_k.*(cos(t/sp.T_max*pi/2-pi/2).^2.*(t<=sp.T_max)+(t>sp.T_max));
case 'AR1'
    % Using first order AR model
    C_max=0.99;
    n_max=sp.Fs*sp.T_max;
    a_T_max=-(1-C_max)^(1/(n_max+1));
    A=filter([1+a_T_max],[1 a_T_max],exp(a_60*t)*sp.A_k,[],1);
otherwise
    error('Bad A_method');
end
% Amplitude coefficient of FM
sp.A_fm=sp.f0*2^(1/12)-sp.f0;
% Frequency coefficient of FM
sp.f_fm=2;
x=A.*cos(2*pi*(sp.f0*t-sp.A_fm/(2*pi*sp.f_fm)*cos(2*pi*sp.f_fm*t+sp.phi_fm))*sp.k_B+sp.phi);
x=sum(x,2);
