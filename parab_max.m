function [x_max,y_max] = parab_max(Px)
% PARAB_MAX - Finds interpolated maximum in Px and returns x and y coordinate.
% This is intended for use on log-spectra to find maxima in between bins.
% x is returned as a fractional index starting at 0 (contrary to common
% MATLAB-style addressing). This is so that it can easily be interpreted as a
% normalized frequency by dividing by the length of Px.
Px=Px(:)';
N=length(Px);
[ma,ind]=max(Px);
indices=mod((ind-1:ind+1)-1,N)+1;
p_ma_est_log_poly=polyfit(indices-1,Px(indices),2);
p_ma_est_log_poly_dx=polyder(p_ma_est_log_poly);
x_max=-p_ma_est_log_poly_dx(2)/(p_ma_est_log_poly_dx(1));
y_max=polyval(p_ma_est_log_poly,x_max);
