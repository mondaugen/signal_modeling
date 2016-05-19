import numpy as np
import matplotlib.pyplot as plt

a0=0
a1=complex('1e-3+0.5j')
N=1024
n=np.arange(N)
j1=complex('1j')
x=np.exp(0+a1*n)
x+=np.random.standard_normal(N)*np.sqrt(0)
plt.figure(1)
plt.plot(n,np.real(x))
w=0.5+0.5*np.cos(2.*np.pi*(np.arange(-N/2.,N/2.))/N)
dw=-np.pi/N*np.sin(2.*np.pi*(np.arange(-N/2.,N/2.))/N)
plt.figure(2)
plt.plot(n,w,n,dw)
# Simply using window
a1_=-np.inner(x,dw)/np.inner(x,w)
print np.imag(a1_)
# Using Fourier transform (bandpass filtering the signal)
#a1__=-np.fft.fft(x*(-w*2.*np.pi*j1*n/N-dw))/np.fft.fft(x*w)
#x_=x[::-1]
#x_=np.r_[x_[-1],x_[:-1]]
Xw=np.fft.fft(x*w)
#Xdw=np.fft.fft(x*dw)
Xdw=np.fft.fft(x*dw)
Xw_=Xw*(-2.*np.pi*j1*n/N)+Xdw
a1__=-Xw_/Xw
print a1__[np.abs(Xw).argmax()]
plt.figure(3)
plt.plot(np.arange(-N/2.,N/2.),np.abs(Xw))
plt.show()
