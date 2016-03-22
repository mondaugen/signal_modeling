function [f,Pf]=spect_max(Px,thresh,F_max,Fs=44100)
% function [f,Pf]=spect_max(Px,thresh,F_max,Fs=44100)
% Assumes Px is full spectrum (redundant for real signals)
N=length(Px);
d=diff(Px);
% only consider maxima greater than threshold
d_max=([0;(d>0)]&[(d<0);0])&(Px>thresh);
idmx=find(d_max);
idmx=idmx(find((idmx/N*Fs)<=F_max));
f=(idmx-1)/N*Fs;
Pf=Px(idmx);
