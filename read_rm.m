function [a,phi,w,psi,mu,t,T,F]=read_rm(path,N,H,Fs)
% [a,phi,w,psi,mu,t,T,F]=READ_RM(L,H,Fs)
%
% path: path to file containing data
% N: length of frame
% H: hopsize
% Fs: sample rate
%
% a: amplitudes
% phi: phases
% w: reassigned frequencies
% psi: estimated frequency slopes
% mu: estimated amplitude slopes
% t: reassigned times
% T: times of each frame (by which to offset t)
% F: frequencies of each frame
f=fopen(path,"r");
if (f == -1)
    error(["Opening file " path]);
end
[data,count]=fread(f,Inf,'double');
fclose(f);
if mod(count,N) != 0,
    error('Length of file must be multiple of N.');
end
N_col=count/N;
data=data(:);
data=reshape(data,[N N_col]);
a=data(:,1:6:end);
phi=data(:,2:6:end);
w=data(:,3:6:end);
psi=data(:,4:6:end);
mu=data(:,5:6:end);
t=data(:,6:6:end);
T=(0:(N_col-1))*H/Fs;
F=(0:(N-1))/N*Fs;
