function [Px] = pgram_hayes(x,N)
%
x=x(:);
if nargin == 1,
    N=1024;
end;
Px=(abs(fft(x,N)).^2)/N;
end;
