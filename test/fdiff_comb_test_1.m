clear;
bw=2;
f=rand(2,1);
f=sort(f)
fdiff=f(2)-f(1);
f_strt=f(1)-floor(f(1)/fdiff)*fdiff;
f_end=f(1)+floor((bw-f(1))/fdiff)*fdiff;
fc=f_strt:fdiff:f_end
