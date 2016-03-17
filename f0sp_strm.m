N_ARGS = 14;
args=argv();
if (nargin != N_ARGS)
    pn=program_name();
    fprintf(stderr,...
    ["usage: " pn " N H W L f0_low f0_high f0_delta K sgm thresh_dB Fs "...
     "F_max P format\n" ...
     "N : FFT size\n" ...
     "H : hop size\n" ...
     "W : window type\n" ...
     "    can be 'hanning', 'hamming'," ...
     " 'blackman', 'bartlett', or 'rectangular'\n" ...
     "L : window length\n" ...
     "f0_low : lowest f0 to test\n" ...
     "f0_high : highest f0 to test\n" ...
     "f0_delta : resolution of f0 in semi-tones\n" ...
     "K : number of harmonics to test for f0 score\n" ...
     "sgm : variance of where harmonics may lie\n" ...
     "thresh_dB: threshold in dB-power below which peaks are not considered\n" ...
     "Fs: samplerate in Hz\n" ...
     "F_max : peaks above this frequency are not considered.\n" ...
     "P : output P pitches with best scores. P < 0, outputs all SCORES.\n" ...
     "format : can be \"raw\" or \"text\". If \"raw\", outputs \"float64\".\n" ...
     "         If \"text\", output human readable values." ...
     "Reads stream of values of type \"float64\" from stdin and ouputs the\n"...
     "  same types to stdout.\n"]);
    exit(-1);
end
N=str2num(args{1});
H=str2num(args{2});
W=args{3};
L=str2num(args{4});
x=zeros(L,1);
f0_low=str2num(args{5});
f0_high=str2num(args{6});
f0_delta=str2num(args{7});
% Fundamentals to test
pitches=(f0_low:f0_delta:f0_high);
f0s=440*2.^((pitches-69)/12);
% Array of scores for each fundamental
K=str2num(args{8});
sgm=str2num(args{9});
thresh_dB=str2num(args{10});
Fs=str2num(args{11});
F_max=str2num(args{12});
P=str2num(args{13});
format=args{14};
fmt=0;
switch format
case "raw"
    fmt=1;
case "text"
    fmt=2;
end
if (!fmt)
    error(sprintf('Bad format %s.',format));
end
% computed this way because we compare with power spectrum
thresh_pow=10^(thresh_dB/10);
% Load in enough samples to fill half of x so that the first frame's centre is
% aligned with the first sample
[tmp,count]=fread(stdin,L/2,'float64');
tmp=tmp(:);
x((L/2+1):(L/2+count))=tmp;
if (count != L/2)
    x((L/2+1+count):L)=zeros(L/2-count,1);
end
do
    s=zeros(length(f0s),1);
    [Px,ENBW]=mper(x,W,N,L);
    d=diff(Px);
    % only consider maxima greater than threshold
    d_max=([0;(d>0)]&[(d<0);0])&(Px>thresh_pow);
    % only consider maxima under nyquist
    d_max=d_max(1:N/2);
    idmx=find(d_max);
    idmx=idmx(find((idmx/N*Fs)<=F_max));
    f=idmx/N*Fs;
    Pf=Px(idmx);
    % compute instantaneous power by dividing out equivalent noise bandwidth of
    % window
    inst_pow=sum(Px)/ENBW;
    for n_=1:length(f0s)
        f0=f0s(n_);
        [s(n_)]=f0_score_puckette(f0,f,Pf,K,sgm);
    end
    if (P>=0)
        [smax,sidx]=sort(s,'descend');
        s=pitches(sidx);
        s=s(1:P);
    end
    switch fmt
    case 1
        % write instantaneous power
        fwrite(stdout,inst_pow,'float64');
        % write f0 scores / pitches
        fwrite(stdout,s,'float64');
    case 2
        % write instantaneous power, then f0 scores / pitches
        outdata=[inst_pow(:);s(:)];
        fprintf(stdout,'%f ',outdata);
        fprintf(stdout,'\n');
    end
    % shift over values
    x(1:(L-H))=x((H+1):L);
    [tmp,count]=fread(stdin,H,'float64');
    tmp=tmp(:);
    x((L-H+1):(L-H+count))=tmp;
until (count != H)
% pad with zeros until buffer full
x((L-H+1+count):L)=zeros(H-count,1);
[Px,ENBW]=mper(x,W,N,L);
d=diff(Px);
% only consider maxima greater than threshold
d_max=([0;(d>0)]&[(d<0);0])&(Px>thresh_pow);
% only consider maxima under nyquist
d_max=d_max(1:N/2);
idmx=find(d_max);
idmx=idmx(find((idmx/N*Fs)<=F_max));
f=idmx/N*Fs;
% compute instantaneous power by dividing out equivalent noise bandwidth of
% window
inst_pow=sum(Px)/ENBW;
for n_=1:length(f0s)
    f0=f0s(n_);
    [s(n_)]=f0_score_puckette(f0,f,Pf,K,sgm);
end
if (P>=0)
    [smax,sidx]=sort(s,'descend');
    s=pitches(sidx);
    s=s(1:P);
end
switch fmt
case 1
    % write instantaneous power
    fwrite(stdout,inst_pow,'float64');
    % write f0 scores / pitches
    fwrite(stdout,s,'float64');
case 2
    % write instantaneous power, then f0 scores / pitches
    outdata=[inst_pow(:);s(:)];
    fprintf(stdout,'%f ',outdata);
    fprintf(stdout,'\n');
end
