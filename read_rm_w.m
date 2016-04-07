function [a,phi,w,psi,mu,t,T,F]=read_rm_w(path,H,Fs)
% [a,phi,w,psi,mu,t,T,F]=READ_RM_W(path,H,Fs)
% WARNING: THIS MIGHT BE BROKEN
%
% path: path to file containing data
% H: hopsize
% Fs: sample rate
%
% Length of each frame read from file
%
% a: amplitudes
% phi: phases
% w: reassigned frequencies
% psi: estimated frequency slopes
% mu: estimated amplitude slopes
% t: reassigned times
% T: times of each frame (by which to offset t)
% F: frequencies of each frame
a=cell();
phi=cell();
w=cell();
psi=cell();
mu=cell();
t=cell();
F=cell();
f=fopen(path,"r");
N_col=6;
if (f == -1)
    error(["Opening file " path]);
end
k=1;
while (~feof(f))
    % Check amount to read
    [M,count]=fread(f,1,'uint64');
    if (count != 1)
        fprintf(stderr,['Error reading file %s\n'...
            'expected %d items, got %d\n'],path,1,count);
        break
    end
    % Read in 6*M worth of doubles
    [data,count]=fread(f,M*N_col,'float64');
    if (count != (M*N_col))
        fprintf(stderr,['Error reading file %s\n'...
            'expected %d items, got %d\n'],path,M*N_col,count);
        break
    end
    data=reshape(data,[M N_col]);
    a{k}=data(:,1);
    phi{k}=data(:,2);
    w{k}=data(:,3);
    psi{k}=data(:,4);
    mu{k}=data(:,5);
    t{k}=data(:,6);
    F{k}=w{k}/(2*pi)*Fs;
    k+=1;
end
T=(0:(k-2))*H/Fs;
fclose(f);
