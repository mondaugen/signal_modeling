function [X,w_,psi_,mu_,t_,X_] = rm_w(x,w,wtype='hanning')
% [X,w_,psi_,mu_,t_,X_] = RM_W(x,w,wtype)
% WARNING: THIS MIGHT BE BROKEN
%
% The reassignment method for computing reassigned frequency, time and frequency
% slope of the STFT of x. From these the frequency modulation and amplitude
% modulation are estimated.
% In this implementation, evaluate the DFT at frequencies w instead of those
% evaluated by the FFT (is therefore a bit slower).
%
% Input arguments:
% x:     a signal of length N. floor(N/2)+1 is the index of x corresponding to time
%        0 (indexing starts at 1).
% w:     The frequencies at which to evaluate, in radians.
% wtype: the window type. Supported types:
%        'hanning'
% 
% Output arguments:
% X:    the original spectrum, divided by sum of window
% w_:   the reassigned frequencies.
% psi_: the estimated frequency modulation parameter (see below).
% mu_:  the esimtated amplitude modulation parameter (see below).
% t_:   the reassigned times.
% X_:   a column vector of length N of complex values representing the initial
%       amplitudes and phases estimated using the reassigned parameters.
%
% The signal x is modeled as a sum of windowed frequency-modulated sinusoids,
% each of the form:
% s(t) = h(t)*a_*exp(mu_*t+j(phi_+w_*t+psi_/2*t^2))
% where the parameters are from above and a_ and phi_ are computed from X_.
x=x(:);
N=length(x);
w=w(:);
n=(1:N)-1;
W=exp(-j*2*pi*w*n);
n_mid=floor(N/2);
t=(1:N)-n_mid-1;
[h,dh,ddh]=rm_h(t,wtype);
ht=h.*t';
dht=dh.*t';
Xht=W*(x.*ht);
Xh=W*(x.*h);
X=Xh/sum(h);
%X=Xh;
Xdh=W*(x.*dh);
Xdht=W*(x.*dht);
Xddh=W*(x.*ddh);
%w=2*pi*(0:(N-1))/N;
%w=w(:);
Xdh_Xh=Xdh./Xh;
w_=w-imag(Xdh_Xh);
t_=real(Xht./Xh);
dw_=imag(Xddh./Xh)-imag((Xdh_Xh).^2);
dt_=real(Xdh.*Xht./(Xh.^2))-real(Xdht./Xh);
psi_=dw_./dt_;
mu_=-real(Xdh_Xh);
dw=-imag(Xdh_Xh);
gam=exp(mu_*t+j*(dw_*t+0.5*psi_*t.^2))*h;
X_=Xh./gam;

function [h,dh,ddh]=rm_h(t,wtype)
% Computes window and window derivatives. Assumes dt=1.
t=t(:);
T=length(t);
switch wtype
    case 'hanning'
        h=cos(pi*t/T).^2;
        dh=-pi/T*sin(2*pi*t/T);
        ddh=-2*(pi/T)^2*cos(2*pi*t/T);
    case 'nutall3'
        h=(0.5-0.5*cos(2*pi*t/T)).^2;
        dh=pi/T*(sin(2*pi*t/T)-0.5*sin(4*pi*t/T));
        ddh=2*(pi/T)^2*(cos(2*pi*t/T)-cos(4*pi*t/T));
    otherwise
        error(sprintf('Bad window: %d',wtype));
end
