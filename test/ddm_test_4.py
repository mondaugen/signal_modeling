import numpy as np
import matplotlib.pyplot as plt

a0=0
a1=complex('1e-3+0.5j')
a2=complex('-1e-6+1e-4j')
N=512
n=np.arange(N)+128
j1=complex('1j')
x=np.exp(0+a1*n+a2*n**2.)
x+=np.random.standard_normal(N)*np.sqrt(0.)
plt.figure(1)
plt.plot(n,np.real(x))
#w=0.5+0.5*np.cos(2.*np.pi*(np.arange(-N/2.,N/2.))/N)
#dw=-np.pi/N*np.sin(2.*np.pi*(np.arange(-N/2.,N/2.))/N)
w=0.5+0.5*np.cos(2.*np.pi*np.arange(-(N-2)/2.,(N-2)/2.+1)/(N-1))
dw=-np.pi/N*np.sin(2.*np.pi*np.arange(-(N-2)/2.,(N-2)/2.+1)/(N-1))
plt.figure(2)
plt.plot(n[:-1],w,n[:-1],dw)
# Simply using window
a1_=-np.inner(x[:-1],dw)/np.inner(x[:-1],w)
print np.imag(a1_)
# Using Fourier transform (bandpass filtering the signal)
#Xw=np.fft.fft(x*w)
#Xp2w=np.fft.fft((x*2.*n)*w)
#Xdw=np.fft.fft(x*dw)
n_=np.arange(N)
Xw=np.fft.fft(x[:-1]*w,N)#*np.exp(-2.*np.pi*1j*n/N)
Xp2w=np.fft.fft((x[:-1]*2.*n[:-1])*w,N)#*np.exp(-2.*np.pi*1j*n/N)
Xdw=np.fft.fft(x[:-1]*dw,N)#*np.exp(-2.*np.pi*1j*n/N)
Xw_=Xw*(-2.*np.pi*j1*n_/N)+Xdw
Xw_A=np.abs(Xw)
Xw_mai0=Xw_A.argmax()
#if (Xw_A[Xw_mai0+1]>Xw_A[Xw_mai0-1]):
#    Xw_mai1=Xw_mai0+1
#else:
#    Xw_mai1=Xw_mai0-1
#A=np.array([
#    [Xw[Xw_mai0],Xp2w[Xw_mai0]],
#    [Xw[Xw_mai1],Xp2w[Xw_mai1]]
#    ])
#b=np.array([
#    [Xw_[Xw_mai0]],
#    [Xw_[Xw_mai1]]
#    ])
Xw_mai1=Xw_mai0+1
Xw_mai_1=Xw_mai0-1
A=np.array([
    [Xw[Xw_mai0],Xp2w[Xw_mai0]],
    [Xw[Xw_mai1],Xp2w[Xw_mai1]],
    [Xw[Xw_mai_1],Xp2w[Xw_mai_1]]
    ])
b=np.array([
    [Xw_[Xw_mai0]],
    [Xw_[Xw_mai1]],
    [Xw_[Xw_mai_1]]
    ])
a=np.linalg.lstsq(A,-b)[0]
print a
plt.figure(3)
plt.plot(np.arange(-N/2.,N/2.),np.abs(Xw))
plt.show()
