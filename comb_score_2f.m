function [fdiff,s]=comb_score_2f(f1,f2,fs,sgm=1,alph=0)
%function [fdiff,s]=comb_score_2f(f1,f2,f,sgm=1,alph=0)
fs=fs(:);
f=[f1;f2];
f=sort(f);
bw=max(f);
fdiff=f(2)-f(1);
f_strt=f(1)-floor(f(1)/fdiff)*fdiff;
f_end=f(1)+floor((bw-f(1))/fdiff)*fdiff;
fc=f_strt:fdiff:f_end;
%s=sum(sum(exp(-(fs-fc).^2/(2*sgm^2)).*((1:length(fc)).^alph)));
s=sum(sum(exp(-(fs-fc).^2/(2*sgm^2)).*fs.^alph));
%s=s/length(fc);
