function [S,T,F,P]=f0sc_read(fpath,N,H,Fs=44100)
%function [S,T,F]=f0sc_read(fpath,N,H,Fs=44100)
% Note that this is reading a cepstrum and so the frequency is calculated as
% Fs/k where k is the bin number (starting at 0)
f=fopen(fpath,'r');
S=fread(f,Inf,'float64');
S=reshape(S,[N+1,length(S)/(N+1)]);
P=S(1,:);
P=P(:);
S=S(2:end,:);
T=(0:(size(S,2)-1))/Fs;
T=T(:);
F=Fs./(0:(size(S,1)-1));
F=F(:);

