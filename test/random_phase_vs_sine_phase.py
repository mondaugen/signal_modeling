import numpy as np
import matplotlib.pyplot as plt
#def D(phi,N):
#    # Dirichlet kernel
#    return np.where(phi==0.,N,np.exp(1j*phi/2.)*np.sin(0.5*N*phi)/np.sin(0.5*phi))
def D(k,N):
    return np.where(k==0,float(N),
            (1.-np.exp(-2.*np.pi*1j*k))/(1.-np.exp(-2.*np.pi*1j*k/float(N))))

#def hann_W(phi,N):
#    # DFT of hann window
#    return (0.5*D(phi,N)
#            +0.25*(D(phi-(2.*np.pi)/N,N)+D(phi+(2.*np.pi)/N,N)))

def hann_W(k,N):
    # DFT of hann window
    return (0.5*D(k,N)-0.25*(D(k-1,N)+D(k+1,N)))

N=16
H=4
k_=3.125
m=np.arange(N)
w=0.5-0.5*np.cos(2.*m/float(N)*np.pi)#np.hanning(N+1)[0:N]
H_max=8
n=np.arange(0,H_max,H)
x_=w*np.exp(1j*2.*np.pi*k_*m/float(N))
print x_
X_=np.fft.fft(x_)
print X_
#phi_x=2*np.pi*m/N
k_x=m
#X_th=hann_W(phi_x-k_/N*2.*np.pi,N)
X_th=hann_W(m-k_,N)
plt.figure(1)
plt.plot(m,np.imag(X_),m,np.imag(X_th),m,np.real(X_),m,np.real(X_th))
x_th=np.fft.ifft(X_th)
plt.figure(2)
plt.plot(m,np.real(x_),m,np.real(x_th))
plt.figure(3)
plt.plot(m,D(m,N),m,np.fft.fft(np.ones(N)))
plt.show()
