# Compare 3 frames 1 bin, 1 frame 3 bins and 3 frames 3 bins
# See how the closeness of two sinusoids disrupts the estimation of amplitude
# slope
# Here we try with a continuous blackman window (see
# test/continuous_blackman_1.py to see how the coefficients were found).
# We also add zero padding

import numpy as np
import sigmod as sm
import matplotlib.pyplot as plt
f0=0.01
f1=0.0155
a0_0=0
a0_1=0#-float('Inf')
a1_0=-1e-6+1j*2.*np.pi*f0
a1_1=-1e-5+1j*2.*np.pi*f1
a2_0=0
a2_1=0
print 'First parameters'
print a0_0
print a1_0
print a2_0
print 'Second parameters'
print a0_1
print a1_1
print a2_1
# Size of non-zero part of window
N_w=4096
# Size of FFT, or window with zero-padding
N=4096
H=1024
L=16000
l=np.arange(L)-H-(N-N_w)/2
x0=np.exp(a0_0+a1_0*l+a2_0*l**2.)
x1=np.exp(a0_1+a1_1*l+a2_1*l**2.)
x=x0+x1
n=np.arange(N)
n_w=np.arange(N_w)
# window coefficients
wc=np.r_[0.358735,0.488305,0.141265,0.011695]
#wc=np.r_[0.5,0.5]
w_=((wc*np.power(-1,np.arange(len(wc))))[:,np.newaxis]
        *np.cos(np.pi*2./N_w*np.outer(np.arange(len(wc)),n_w)))
w=np.sum(w_,0)
dw_=((2.*np.pi/N_w*wc[1:]*np.arange(1,len(wc))
        *np.power(-1,1+np.arange(1,len(wc))))[:,np.newaxis]
    *np.sin(np.pi*2./N_w*np.outer(np.arange(1,len(wc)),n_w)))
dw=np.sum(dw_,0)
# zero-padded windows
w_p=np.r_[np.zeros((N-N_w)/2),w,np.zeros((N-N_w)/2)]
dw_p=np.r_[np.zeros((N-N_w)/2),dw,np.zeros((N-N_w)/2)]
#w=0.5-0.5*np.cos(2.*np.pi*n/N)
#dw=np.pi/N*np.sin(2.*np.pi*n/N)
a=sm.ddm_p2_3_1(x,H,w_p,dw_p,4)
print a
a_=sm.ddm_p2_1_3(x[H:],w_p,dw_p)
print a_
a__=sm.ddm_p2_3_3(x,H,w_p,dw_p,4)
print a__

X=np.fft.fft(x[H:H+N]*w_p)
plt.plot(n,20.*np.log10(np.abs(X)/np.sum(w_p)))
plt.show()
