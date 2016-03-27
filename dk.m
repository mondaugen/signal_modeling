function [y]=dk(theta,N)
% [y]=dk(theta,N)
% Compute the dirichlet kernel of order N.
% theta is in radians
% Multiplied by complex exponential to centre around N/2 in time-domain
y=exp(-j*N*theta/2).*sin((N/2+1/2)*theta)./sin(theta/2);
% Remove NaNs because function is defined for theta=0
y(find(theta==0))=N;
