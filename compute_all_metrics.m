function [Pxm,PX,ENBW_PX]=compute_all_metrics(opt)
% [Pxm]=compute_all_metrics(fname,opt)
%
% Load in file of float64 and compute a set of metrics on each frame.
% opt is a set of options where among other things, you specify the file.
%
% Pxm : A cell containing as many entries as there are columns in PX, each
%       containing a structure of information about that frame.
% PX : A matrix containing the powespectrum of each frame in its columns
% ENBW_PX : The equivalent noise bandwidth of the window used to compute the
%           power spectrum, useful for computing the total power in a frame.

if (nargin < 1)
    opt=struct();
end

% Default fields and values
def_flds=struct(
    'fname','/tmp/test.f64',
    'Fs',16000,
% Window size, for buffering
    'L',1024,
% Hop size, for buffering
    'H',256,
% DFT size
    'N',2048,
% Window to use when computing power spectrum
    'win','hanning',
% Number of frames to compute
    'K',20,
% Peak picking algorithm parameters
    'ppvp_a',0.1,
    'ppvp_c',0.5,
    'ppvp_q',2,
% Window to use with reassignment method.
    'rm_win','blackman',
% f0s to compute scores for, specified in MIDI pitch (69 == A440)
    'f0s_pchs',(36:84)',
% generates harmonics no higher than this frequency
    'f0s_max_f',5000,
% harmonic comb tooth spread
    'f0s_sd_f',20,
% rate of harmonic decay db/octave
    'f0s_h_slp',-6
);

% Check what fields are present and replace with defaults
for fn=fieldnames(def_flds)'
    fn=char(fn);
    if(~isfield(opt,fn))
        opt.(fn)=def_flds.(fn);
    end
end

fi=fopen(opt.fname,'r');
x=fread(fi,Inf,'float64');
x+=randn(length(x),1)*0.01;
fclose(fi);
% Put file into frames
X=buffer_ne(x,opt.L,opt.H);
% Compute powerspectrum on each frame
[PX,ENBW_PX]=mper(X,opt.win,opt.N,opt.L-1);
% Cell to store information
Pxm=cell();
for n=1:opt.K
    % Struct to store values
    st=struct();
    Xn=PX(:,n);
    % Find extrema in power spectrum
    [st.Xm,st.Xmi]=lextrem(Xn,'max');
    % Interpolate to obtain better values 
    % (only works when DFT size at least twice window length)
    [st.Xm_,st.Xmi_]=lmax_parab(p2db(Xn),st.Xmi);
    % Filter out erroneous maxima using algorithm
    % (the resulting indices can be used to look up the maxima's values in
    % st.Xmi_)
    [st.ex_i]=ppvp(st.Xm_,opt.ppvp_a,'max',opt.ppvp_c,opt.ppvp_q);
    % Frequencies at estimated maxima
    st.fm_=(st.Xmi_(st.ex_i)-1)/length(Xn)*opt.Fs;
    % Compute reassigned parameters
    [st.X_r,st.w_r,st.psi_r,st.mu_r,st.t_r,st.X_r_]=rm(X(:,n),opt.rm_win);
    % RM computes power spectrum differently so we compute its maxima
    % Find extrema in spectrum, make power spectrum
    [st.X_r_m,st.X_r_mi]=lextrem(abs(st.X_r).^2,'max');
    % Filter out erroneous maxima using algorithm
    % (the resulting indices can be used to look up the maxima's values in
    % st.X_r_m)
    [st.X_r_ex_i]=ppvp(st.X_r_m,opt.ppvp_a,'max',opt.ppvp_c,opt.ppvp_q);
    % Convert pitches to Hz
    opt.f0s=midi2hz(opt.f0s_pchs);
    st.f0_s=zeros(size(opt.f0s,1),length(st.X_r_ex_i));
    % Computes f0 scores using RM maxima
    p=st.w_r(st.X_r_mi(st.X_r_ex_i))/(2*pi)*opt.Fs;
    % store reassigned frequencies
    st.f_r=p;
    p=p(:)';
    % (TODO: Probably better would be to use a value interpolated for the
    % reassigned frequency)
    Ap=abs(st.X_r)(st.X_r_mi(st.X_r_ex_i));
    % Store (interpolated) frequencies
    st.A_r=Ap;
    for _k=1:size(opt.f0s,1)
        st.f0_s(_k,:)=p_f0_score(opt.f0s(_k),p,Ap,opt.f0s_max_f,opt.f0s_sd_f,opt.f0s_h_slp);
    end
    Pxm{n}=st;
end
