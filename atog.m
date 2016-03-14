function [gamma]=atog(a)
% Find the reflection coefficients gamma from the polynomial coefficients a
a=a(:);
p=length(a);
a=a(2:p)/a(1);
gamma(p-1)=a(p-1);
for j=p-1:-1:2,
    a=(a(1:j-1)-gamma(j)*flipud(conj(a(1:j-1))))./ ...
        (1-abs(gamma(j))^2);
    gamma(j-1)=a(j-1);
end;
gamma=gamma(:);
