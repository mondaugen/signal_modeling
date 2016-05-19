import numpy as np
import matplotlib.pyplot as plt
import sigmod as sm

a0=0
a1=complex('1e-3+0.5j')
a2=complex('-1e-6+1e-4j')
N_w=511
H=128
N=2*H+N_w+1
n=np.arange(N)-H
j1=complex('1j')
x=np.exp(0+a1*n+a2*n**2.)
x+=np.random.standard_normal(N)*np.sqrt(0)
#plt.figure(1)
#plt.plot(n,np.real(x))
w=0.5+0.5*np.cos(2.*np.pi*(np.arange(-(N_w-1)/2.,(N_w-1)/2.+1))/N_w)
dw=-np.pi/N*np.sin(2.*np.pi*(np.arange(-(N_w-1)/2.,(N_w-1)/2.+1))/N_w)
a_9=sm.ddm_p2_3_3(x,H,w,dw,3)
a_3=sm.ddm_p2_1_3(x,w,dw)
print 'True:'
print np.c_[np.r_[a1,a2]]
print 'Estimated with 9 atoms:'
print a_9
print 'Estimated with 3 atoms:'
print a_3
#plt.figure(2)
#plt.plot(n,w,n,dw)
# Simply using window
# Using Fourier transform (bandpass filtering the signal)
#plt.figure(3)
#plt.plot(np.arange(-N/2.,N/2.),np.abs(Xw))
#plt.show()

