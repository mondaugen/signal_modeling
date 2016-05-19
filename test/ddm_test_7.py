import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

H=128
a0=1j
a1=complex('1e-3+0.5j')
a2=complex('-1e-5+1e-4j')
N_w=512
n=np.arange(N_w+2*H)-H
j1=complex('1j')
x=np.exp(a0+a1*n+a2*n**2.)
x+=np.random.standard_normal(N_w+2*H)*np.sqrt(0.1)
w=0.5+0.5*np.cos(2.*np.pi*(np.arange(-N_w/2.,N_w/2.)/N_w))
dw=-np.pi/N_w*np.sin(2.*np.pi*(np.arange(-N_w/2.,N_w/2.)/N_w))
a_9=sm.ddm_p2_3_3(x,H,w,dw,2)
print 'True:'
print np.c_[np.r_[a0,a1,a2]]
print 'Estimated with 9 atoms:'
print a_9
#plt.figure(2)
#plt.plot(n,w,n,dw)
# Simply using window
# Using Fourier transform (bandpass filtering the signal)
#plt.figure(3)
#plt.plot(np.arange(-N/2.,N/2.),np.abs(Xw))
#plt.show()

