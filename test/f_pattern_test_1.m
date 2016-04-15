% Take a sequence of harmonics, round, remove duplicate elements, sort 
% Divide all harmonics by the number that was rounded to (e.g., round to 0.05->
% to round 1.2345 we do round(1.2345/0.05)*0.05).
% This makes a "signature". Using this signature, can we look up a
% plausible combination of fundamental frequencies giving this sequence of
% harmonics?
clear;
p_max=12;
p=(0:p_max)';
k=1:10;
A=2.^(p/12)*k;
% Round
A=round(A/0.05)*0.05;
R=cell();
S=cell();
for n=2:length(p)
    m=n-1;
    S{m}=unique([A(1,:),A(m,:)]);
    R{m}=S{m}/0.05;
end
