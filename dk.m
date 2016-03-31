function [y]=dk(theta,N)
% [y]=dk(theta,N)
% Compute the dirichlet kernel of order N.
% theta is in radians
y=sin((N/2+1/2)*theta)./sin(theta/2);
% Remove NaNs because function is defined for theta=0
% Some care should be taken in the future to handle values of theta other than 0
% that give small values of sin but so far no real problems have been observed.
y(find(theta==0))=(N+1);
if (mod(N,2)==0)
    % N even
    idx=find(mod(theta,2*pi)==0);
    y(idx)=(N+1)*(-1).^(theta(idx)/(2*pi));
elseif (mod(N,2)==1)
    % N odd
    y(find(mod(theta,2*pi)==0))=(N+1);
end
