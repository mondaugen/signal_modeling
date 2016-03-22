% name of soundfile and log file
filepre='/tmp/2016-03-22T18:47:12,503968000-0400_mix';
sndfile_p=[filepre '.f64le'];
logfile_p=[filepre '.INFO'];
%sndfile_f=fopen(sndfile_p,'r');
logfile_f=fopen(logfile_p,'r');
logfile_t=fread(logfile_f);
fclose(logfile_f);
logfile_t=char(logfile_t(:)');
% extra newline to help parsing
logfile_t=["\n" logfile_t];
pitches=eval(['['...
    strsplit(logfile_t(regexp(logfile_t,"\npitches=")+9:end)){1}...
    ']']);
files=eval(['['...
    strsplit(logfile_t(regexp(logfile_t,"\nfiles=")+7:end)){1}...
    ']']);
f0sc_strm_f(8192,1024,'blackman',4096,24,96,0.5,44100,-1,'raw',sndfile_p);
[S,T,F]=f0sc_read([sndfile_p '.f0sc'],8192,1024,44100);
imagesc(T,(0:size(S,1))/44100,S);
%imagesc(S);
% function f0sc_strm_f(N,H,W,L,pch_low,pch_high,c,Fs,P,format,datapath)
%logfile_t=strsplit(logfile_t,"\n");
