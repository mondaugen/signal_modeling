N_ARGS=13;
argv=argv();
if length(argv) != N_ARGS
    fprintf(stderr,
    ['Arguments are:\n'...
    'sample-rate\n'...
    'length-samples\n'...
    'fundamental-Hz\n'...
    'number-harmonics\n'...
    'inharmonicity-coefficient\n'...
    '60dB-harmonic-number\n'...
    'initial-phase\n'...
    'initial-fm-phase\n'...
    't-60dB\n'...
    'attack-max-time\n'...
    'attack-method {FOF,AR1}\n'...
    'amplitude-fm\n'...
    'frequency-fm\n']);
    exit(-1);
end

opt=struct();
% Sampling rate in Hz
opt.Fs=str2num(argv{1});
% Length of signal
opt.N=str2num(argv{2});
% Fundamental
opt.f0=str2num(argv{3});
% Number of harmonics
opt.K=str2num(argv{4});
% Inharmonicity coefficient
%opt.B=exp(2.54*log(opt.f0)-24.6); % For now, the just noticeable inharmonic coefficient
opt.B=str2num(argv{5});
% Amplitude based on harmonic number
opt.A_k_60=str2num(argv{6}); % The harmonic number that is 60dB lower than the first
% Initial phase
opt.phi=str2num(argv{7});
% Initial FM phase
opt.phi_fm=str2num(argv{8});
% Time in seconds until amplitude of partial has dropped by 60dB
%opt.T_60=1;
opt.T_60=str2num(argv{9});
% Time in seconds to when attack function reaches maximum
opt.T_max=str2num(argv{10});
opt.A_method=argv{11};
opt.A_fm=str2num(argv{12});
opt.f_fm=str2num(argv{13});
x=harm_sines(opt);
fwrite(stdout,x,'float64');
