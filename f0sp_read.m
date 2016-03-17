function [S,P,pitches]=f0sp_read(path,f0_low,f0_high,f0_delta)
% [S,P]=f0sp_read(path,f0_low,f0_high,f0_delta)
% Read a file of f0 scores.
%
% S : The scores for each f0.
% P : The signal power for each frame
% path : path to the scores
% f0_low : the lowest pitch considered. 69 = 440Hz
% f0_high : the highest pitch considered
% f0_delta : the resolution of pitches in semitones
f=fopen(path,'r');
data=fread(f,Inf,'float64');
fclose(f);
pitches=f0_low:f0_delta:f0_high;
data=reshape(data,[length(pitches)+1,length(data)/(length(pitches)+1)]);
P=data(1,:);
S=data(2:end,:);
