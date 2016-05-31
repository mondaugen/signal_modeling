# Compare 3 frames 1 bin, 1 frame 3 bins and 3 frames 3 bins
# See how the closeness of two sinusoids disrupts the estimation of amplitude
# slope
import numpy as np
import sigmod as sm
import matplotlib.pyplot as plt
f0=0.01
f1=0.02
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
N=512
H=128
L=2048
l=np.arange(L)-H
x0=np.exp(a0_0+a1_0*l+a2_0*l**2.)
x1=np.exp(a0_1+a1_1*l+a2_1*l**2.)
x=x0+x1
n=np.arange(N)
w=0.5-0.5*np.cos(2.*np.pi*n/N)
dw=np.pi/N*np.sin(2.*np.pi*n/N)
a=sm.ddm_p2_3_1(x,H,w,dw,2)
print a
a_=sm.ddm_p2_1_3(x[H:],w,dw)
print a_
a__=sm.ddm_p2_3_3(x,H,w,dw,2)
print a__

X=np.fft.fft(x[H:H+N]*w)
plt.plot(n,20.*np.log10(np.abs(X)/np.sum(w)))
plt.show()
