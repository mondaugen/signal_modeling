function [S,T,F] = read_half_spect(path,N,H,Fs)
% [S,T,F] = READ_HALF_SPECT(path,N,H,Fs)
%
% path: path to spectrogram stored as array of "float64"
% N:    number of bins per frame
% H:    hop size in samples
% Fs:   sample rate
%
% S will be a matrix of N/2 rows and length(file)/H columns
m=1;
f=fopen(path,"r");
if (f == -1)
    error(["Opening file " path]);
end
[data,count]=fread(f,Inf,'double');
fclose(f);
% if count not multiple of N, ignore remainder
if mod(count,N) != 0,
    error('Length of file must be multiple of N.');
end
N_col=(count/2)/N;
data=data(:);
data_real=data(1:2:end);
data_imag=data(2:2:end);
data_real=reshape(data_real,[N,N_col]);
data_imag=reshape(data_imag,[N,N_col]);
S=data_real+1i.*data_imag;
S=S(1:(N/2),:);
T=(0:(N_col-1))*H/Fs;
F=(0:(N/2-1))/N*Fs;
