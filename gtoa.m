function [a]=gtoa(gamma)
% gtoa - Find the polynomial coefficients a from the reflection coefficients
% gamma.
a=1;
gamma=gamma(:);
p=length(gamma);
for j=2:p+1,
    a=[a;0]+gamma(j-1)*[0;conj(flipud(a))];
end;
