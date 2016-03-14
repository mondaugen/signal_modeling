function [A_]=ar_stabilizer(A),
% Stabilize A (make it minimum phase) by finding its roots and then replacing
% those that are outside the unit circle with their conjugated inverse. The
% polynomial is then put back together. A and A_ are filter coefficients in the
% form that is accepted for the A argument of the "filter" function.
A=A(:);
R=roots(A);
R(abs(R)>1)=conj(1./R(abs(R)>1));
A_=[1];
for n=(1:length(R)),
    A_=conv(A_,[1 -R(n)]);
end;
A_=A_(:);
