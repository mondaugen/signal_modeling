function [X,w_,psi_,mu_,t_,X_] = rm(x,wtype='hanning')
% [X_,w_,psi_,mu_,t_] = RM(x,wtype)
%
% The reassignment method for computing reassigned frequency, time and frequency
% slope of the STFT of x. From these the frequency modulation and amplitude
% modulation are estimated.
%
% Input arguments:
% x:     a signal of length N. floor(N/2)+1 is the index of x corresponding to time
%        0 (indexing starts at 1).
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
n_mid=floor(N/2);
t=(1:N)-n_mid-1;
[h,dh,ddh]=rm_h(t,wtype);
ht=h.*t';
dht=dh.*t';
Xht=fft(x.*ht);
Xh=fft(x.*h);
X=Xh/sum(h);
Xdh=fft(x.*dh);
Xdht=fft(x.*dht);
Xddh=fft(x.*ddh);
w=2*pi*(0:(N-1))/N;
w=w(:);
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
    case {'blackman3-min','blackman'}
        a=[0.42323,0.49755,0.07922];
        if (mod(T,2)==0)
            % T is even
            D=T;
        else
            D=T-1;
        end
        M=length(a);
        m=0:(M-1);
        h=sum(a.*cos(2*pi/D*t*m),2);
        dh=-1*sum((2*pi/D*m(2:end)).*a(2:end).*sin(2*pi/D*t*m(2:end)),2);
        ddh=-1*sum((2*pi/D*m(2:end)).^2.*a(2:end).*cos(2*pi/D*t*m(2:end)),2);
    case 'nutall3'
        h=(0.5-0.5*cos(2*pi*t/T)).^2;
        dh=pi/T*(sin(2*pi*t/T)-0.5*sin(4*pi*t/T));
        ddh=2*(pi/T)^2*(cos(2*pi*t/T)-cos(4*pi*t/T));
    otherwise
        error(sprintf('Bad window: %d',wtype));
end
